#!/usr/bin/env python3
"""
joint_fp_engine.py - Production Inference Engine for Joint 2215-bit Fingerprint Model:
- Pure PyTorch Transformer Model (MS2 -> Joint 2215-bit ECFP2048+MACCS167)
- Downstream QSAR environmental toxicity scoring (Random Forest on 2215-d vectors)
- Decoupled functional group parsing (MACCS key annotations)
- Self-contained CLI and importable Python API for NTSminer R wrapper integration.
"""

import sys
if sys.version_info < (3, 8):
    sys.exit("Error: Python 3.8 or higher is required for IDSL_MINT joint fingerprint engine.")

import os
os.environ["PYTORCH_CUDA_ALLOC_CONF"] = "expandable_segments:True"

import json
import argparse
from pathlib import Path
from typing import List, Dict, Any, Tuple, Optional

import torch
import torch.nn as nn
import numpy as np

# --------------------------------------------------------------------------
# 1. Self-contained MS2FP Neural Architecture
# --------------------------------------------------------------------------
class MS2FP_Model(nn.Module):
    def __init__(self,
                 NUM_MZ_VOCABS: int = 148003,
                 D_MODEL: int = 256,
                 MZ_EMBED_MAX_NORM: int = 2,
                 MZ_DROPOUT: float = 0.0,
                 MAX_NUMBER_BITS: int = 280,
                 NUM_FINGERPRINT_VOCABS: int = 2218,
                 NUM_HEAD: int = 8,
                 NUM_ENCODER_LAYER: int = 4,
                 NUM_DECODER_LAYER: int = 4,
                 TRANSFORMER_DROPOUT: float = 0.0,
                 TRANSFORMER_ACTIVATION: str = 'relu'):
        super().__init__()
        self.D_model = torch.tensor(D_MODEL)
        self.max_number_bits = MAX_NUMBER_BITS
        self.vocab_size = NUM_FINGERPRINT_VOCABS

        self.mzTokenEmbedding = nn.Embedding(
            num_embeddings = NUM_MZ_VOCABS + 1,
            embedding_dim = D_MODEL,
            max_norm = MZ_EMBED_MAX_NORM
        )
        self.mzDropout = nn.Dropout(p = MZ_DROPOUT)
        
        self.nBitsEmbedding = nn.Embedding(
            num_embeddings = NUM_FINGERPRINT_VOCABS + 1,
            embedding_dim = D_MODEL,
            padding_idx = 0,
            max_norm = 2
        )

        self.mzTransformer = nn.Transformer(
            d_model = D_MODEL,
            nhead = NUM_HEAD,
            num_encoder_layers = NUM_ENCODER_LAYER,
            num_decoder_layers = NUM_DECODER_LAYER,
            dim_feedforward = 4 * D_MODEL,
            dropout = TRANSFORMER_DROPOUT,
            activation = TRANSFORMER_ACTIVATION,
            batch_first = False
        )

        self.Perception = nn.Linear(in_features = D_MODEL, out_features = NUM_FINGERPRINT_VOCABS, bias = True)

    def PeakEmbedding(self, MZ_Tokens: torch.Tensor, INT: torch.Tensor) -> torch.Tensor:
        z = self.mzTokenEmbedding(MZ_Tokens) * torch.sqrt(self.D_model.to(MZ_Tokens.device))
        device = z.device
        idx2i = torch.arange(start = 0, end = self.D_model, step = 2, device = device)
        idx1i = torch.arange(start = 1, end = self.D_model, step = 2, device = device)

        zPosEnc = torch.zeros_like(z, device = device)
        zPosEnc[:, :, idx2i] = torch.sin(idx2i * INT / (10000**(2 * idx2i / self.D_model)))
        zPosEnc[:, :, idx1i] = torch.cos(idx1i * INT / (10000**(2 * idx1i / self.D_model)))
        return self.mzDropout(z + zPosEnc)

    def forward(self, MZ_Tokens: torch.Tensor, INT: torch.Tensor, FP_Tokens: torch.Tensor, tgt_key_padding_mask = None) -> torch.Tensor:
        Source = self.PeakEmbedding(MZ_Tokens, INT).permute(1, 0, 2)
        Target = self.nBitsEmbedding(FP_Tokens).permute(1, 0, 2)

        T = Target.size(0)
        tgt_mask = torch.triu(torch.full((T, T), float('-inf'), device = Target.device), diagonal = 1)

        TransformedMZ = self.mzTransformer(
            src = Source,
            tgt = Target,
            src_mask = None,
            tgt_mask = tgt_mask,
            memory_mask = None,
            src_key_padding_mask = None,
            tgt_key_padding_mask = tgt_key_padding_mask,
            memory_key_padding_mask = None
        )
        return self.Perception(TransformedMZ.permute(1, 0, 2))

    def encode_spectrum(self, MZ_Tokens: torch.Tensor, INT: torch.Tensor) -> torch.Tensor:
        """Encodes MS2 spectrum into 256-dimensional semantic embedding via Transformer Encoder."""
        Source = self.PeakEmbedding(MZ_Tokens, INT).permute(1, 0, 2)
        memory = self.mzTransformer.encoder(Source)
        embedding = memory.mean(dim=0)
        return nn.functional.normalize(embedding, p=2, dim=-1)


# --------------------------------------------------------------------------
# 2. Key MACCS Feature Interpretation Dictionary
# --------------------------------------------------------------------------
MACCS_DESCRIPTIONS = {
    155: "Fluorine atom (F)",
    146: "Aliphatic C-F bond",
    140: "Carboxylic acid (-COOH) / Carbonyl",
    131: "Sulfonate (-SO3) / Sulfate",
    120: "Ether linkage (C-O-C)",
    143: "Aliphatic chain (C-C-C)",
    164: "Saturated carbon chain",
    139: "Halogenated aliphatic moiety",
    124: "Oxygen-containing polar headgroup",
    115: "Nitrogen heteroatom / Amide",
    105: "Aromatic / Heteroaromatic ring",
    80:  "Sulfur heteroatom"
}

# --------------------------------------------------------------------------
# 3. Spectrum Preprocessing & Tokenization
# --------------------------------------------------------------------------
def tokenize_spectrum(precursor_mz: float,
                      peaks: List[Tuple[float, float]],
                      min_mz: float = 20.0,
                      max_mz: float = 1500.0,
                      interval_mz: float = 0.01,
                      max_peaks: int = 512) -> Tuple[np.ndarray, np.ndarray]:
    """Tokenizes precursor and MS2 peak list into discrete index bins and intensities."""
    # Precursor token
    try:
        precursor_mz = float(precursor_mz)
        if not np.isfinite(precursor_mz):
            precursor_mz = min_mz
    except (ValueError, TypeError):
        precursor_mz = min_mz

    if precursor_mz < min_mz or precursor_mz > max_mz:
        p_token = int(round((max(min_mz, min(max_mz, precursor_mz)) - min_mz) / interval_mz)) + 1
    else:
        p_token = int(round((precursor_mz - min_mz) / interval_mz)) + 1

    valid_peaks = []
    for p in (peaks or []):
        if not isinstance(p, (list, tuple)) or len(p) < 2:
            continue
        mz, intensity = p[0], p[1]
        if mz is None or intensity is None:
            continue
        try:
            mz_f = float(mz)
            int_f = float(intensity)
            if min_mz <= mz_f <= max_mz and np.isfinite(mz_f) and np.isfinite(int_f) and int_f > 0:
                tok = int(round((mz_f - min_mz) / interval_mz)) + 1
                valid_peaks.append((tok, int_f))
        except (ValueError, TypeError):
            continue

    # Sort peaks by intensity descending and truncate to max_peaks
    valid_peaks.sort(key=lambda x: x[1], reverse=True)
    valid_peaks = valid_peaks[:max_peaks]

    # Normalize intensity to 0..100
    if valid_peaks:
        max_int = max(p[1] for p in valid_peaks)
        norm_factor = 100.0 / max(max_int, 1e-6)
    else:
        norm_factor = 1.0

    mz_tokens = [p_token] + [p[0] for p in valid_peaks]
    int_tokens = [100.0] + [min(100.0, p[1] * norm_factor) for p in valid_peaks]

    return np.array(mz_tokens, dtype=np.int32), np.array(int_tokens, dtype=np.float32)


# --------------------------------------------------------------------------
# 4. Joint Fingerprint & Toxicity Inference Pipeline
# --------------------------------------------------------------------------
class JointFPInferencePipeline:
    def __init__(self,
                 weights_path: str,
                 config_path: str,
                 tox_model_path: Optional[str] = None,
                 device: str = "auto"):
        self.device = torch.device("cuda" if (device == "cuda" or (device == "auto" and torch.cuda.is_available())) else "cpu")
        
        # Load configuration
        with open(config_path, "r", encoding="utf-8") as f:
            self.config = json.load(f)
            
        arch = self.config["architecture"]
        self.model = MS2FP_Model(
            NUM_MZ_VOCABS = arch["num_mz_vocabs"],
            D_MODEL = arch["d_model"],
            MAX_NUMBER_BITS = arch["max_sequence_length"],
            NUM_FINGERPRINT_VOCABS = arch["vocab_size"],
            NUM_HEAD = arch["num_heads"],
            NUM_ENCODER_LAYER = arch["num_encoder_layers"],
            NUM_DECODER_LAYER = arch["num_decoder_layers"]
        ).to(self.device)

        # Load weights
        state_dict = torch.load(weights_path, map_location=self.device, weights_only=False)
        self.model.load_state_dict(state_dict)
        self.model.eval()

        # Load Downstream QSAR Toxicity Model
        self.tox_clf = None
        if tox_model_path and Path(tox_model_path).exists():
            try:
                import joblib
                tox_bundle = joblib.load(tox_model_path)
                self.tox_clf = tox_bundle.get("model", tox_bundle) if isinstance(tox_bundle, dict) else tox_bundle
            except Exception as e:
                self.tox_clf = None

    def predict_spectrum(self, precursor_mz: float, peaks: List[Tuple[float, float]], max_len: int = 280) -> Dict[str, Any]:
        """Runs autoregressive greedy inference for a single MS2 spectrum."""
        arch = self.config["architecture"]
        mz_toks, int_toks = tokenize_spectrum(
            precursor_mz, peaks,
            min_mz = arch["min_mz"],
            max_mz = arch["max_mz"],
            interval_mz = arch["interval_mz"],
            max_peaks = arch["max_peaks"]
        )

        t_mz = torch.tensor(mz_toks, dtype=torch.int, device=self.device).unsqueeze(0)
        t_int = torch.tensor(int_toks, dtype=torch.float32, device=self.device).unsqueeze(0).unsqueeze(2)

        # Greedy decoding starting with [SOS] = 1
        curr_tokens = torch.tensor([[1]], dtype=torch.long, device=self.device)
        
        with torch.no_grad():
            for _ in range(max_len):
                logits = self.model(t_mz, t_int, curr_tokens)
                next_tok = logits[:, -1, :].argmax(dim=-1, keepdim=True)
                tok_val = next_tok.item()
                if tok_val == 2:  # [EOS]
                    break
                curr_tokens = torch.cat([curr_tokens, next_tok], dim=1)

        predicted_tokens = curr_tokens.squeeze(0).cpu().numpy()
        
        # Extract active on-bits (tokens >= 3 map to bits 0..2214)
        active_bits = sorted(list(set(int(t - 3) for t in predicted_tokens if 3 <= t < 2218)))
        
        # Split into ECFP2048 (bits 0..2047) and MACCS167 (bits 2048..2214)
        ecfp_bits = [b for b in active_bits if b < 2048]
        maccs_bits = [b - 2048 for b in active_bits if b >= 2048]

        # Construct binary 2215-bit vector
        fp_vec = np.zeros(2215, dtype=np.uint8)
        fp_vec[active_bits] = 1

        # Downstream Toxicity Scoring
        tox_prob = None
        tox_risk = "Unknown"
        if self.tox_clf is not None:
            try:
                feat_matrix = fp_vec.reshape(1, -1)
                if hasattr(self.tox_clf, "predict_proba"):
                    tox_prob = float(self.tox_clf.predict_proba(feat_matrix)[0, 1])
                else:
                    tox_prob = float(self.tox_clf.predict(feat_matrix)[0])
                    
                if tox_prob >= 0.70:
                    tox_risk = "High"
                elif tox_prob >= 0.40:
                    tox_risk = "Moderate"
                else:
                    tox_risk = "Low"
            except Exception:
                tox_prob = None
                tox_risk = "Unknown"

        # Functional group interpretations
        detected_groups = [desc for bit, desc in MACCS_DESCRIPTIONS.items() if bit in maccs_bits]

        return {
            "precursor_mz": float(precursor_mz),
            "n_peaks": len(peaks),
            "n_total_bits": len(active_bits),
            "n_ecfp_bits": len(ecfp_bits),
            "n_maccs_bits": len(maccs_bits),
            "joint_fp_on_bits": active_bits,
            "ecfp_on_bits": ecfp_bits,
            "maccs_on_bits": maccs_bits,
            "toxicity_score": round(tox_prob, 4) if tox_prob is not None else None,
            "toxicity_risk": tox_risk,
            "detected_functional_groups": detected_groups
        }

    def encode_spectrum(self, precursor_mz: float, peaks: List[Tuple[float, float]]) -> List[float]:
        """Extracts 256-dimensional L2-normalized semantic embedding via Transformer Encoder."""
        arch = self.config["architecture"]
        mz_toks, int_toks = tokenize_spectrum(
            precursor_mz, peaks,
            min_mz = arch["min_mz"],
            max_mz = arch["max_mz"],
            interval_mz = arch["interval_mz"],
            max_peaks = arch["max_peaks"]
        )

        t_mz = torch.tensor(mz_toks, dtype=torch.int, device=self.device).unsqueeze(0)
        t_int = torch.tensor(int_toks, dtype=torch.float32, device=self.device).unsqueeze(0).unsqueeze(2)

        with torch.no_grad():
            emb = self.model.encode_spectrum(t_mz, t_int)

        return [round(float(v), 6) for v in emb.squeeze(0).cpu().numpy().tolist()]

    def encode_batch(self, spectra_list: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
        results = []
        for item in spectra_list:
            feat_id = item.get("id") or item.get("name") or ("FT_%04d" % (len(results) + 1))
            precursor_mz = item.get("precursor_mz") or item.get("mz") or 0.0
            peaks = item.get("peaks") or []
            emb = self.encode_spectrum(precursor_mz, peaks)
            results.append({
                "id": str(feat_id),
                "precursor_mz": float(precursor_mz),
                "n_peaks": len(peaks),
                "embedding": emb
            })
        return results

    def predict_batch(self, spectra_list: List[Dict[str, Any]], include_embedding: bool = False) -> List[Dict[str, Any]]:
        results = []
        for item in spectra_list:
            feat_id = item.get("id") or item.get("name") or ("FT_%04d" % (len(results) + 1))
            precursor_mz = item.get("precursor_mz") or item.get("mz") or 0.0
            peaks = item.get("peaks") or []
            pred = self.predict_spectrum(precursor_mz, peaks)
            pred["id"] = str(feat_id)
            if include_embedding:
                pred["embedding"] = self.encode_spectrum(precursor_mz, peaks)
            results.append(pred)
        return results


# --------------------------------------------------------------------------
# 5. SMILES 2215-bit Fingerprint Calculator (RDKit Backend)
# --------------------------------------------------------------------------
def compute_smiles_joint_fp(smiles: str) -> Optional[List[int]]:
    """Computes 2215-bit joint fingerprint active on-bits (ECFP2048 + MACCS167 offset 2048)."""
    try:
        from rdkit import Chem
        from rdkit.Chem import AllChem, MACCSkeys
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return None
        fp_ecfp = AllChem.GetMorganFingerprintAsBitVect(mol, radius=2, nBits=2048)
        ecfp_bits = list(fp_ecfp.GetOnBits())
        fp_maccs = MACCSkeys.GenMACCSKeys(mol)
        maccs_bits = [int(b + 2048) for b in fp_maccs.GetOnBits()]
        return sorted(list(set(ecfp_bits + maccs_bits)))
    except Exception:
        return None


def compute_smiles_batch(items: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
    results = []
    for it in items:
        cid = it.get("id", it.get("name", "CAND_%d" % (len(results) + 1)))
        smi = it.get("smiles", it.get("SMILES", ""))
        bits = compute_smiles_joint_fp(smi) if smi else None
        results.append({
            "id": cid,
            "smiles": smi,
            "joint_fp_on_bits": bits,
            "n_bits": len(bits) if bits is not None else 0
        })
    return results


def extract_murcko_scaffold(smiles: str) -> Dict[str, Any]:
    """Extracts Bemis-Murcko scaffold SMILES and generic framework using RDKit."""
    if not smiles or not isinstance(smiles, str) or not smiles.strip():
        return {"scaffold_smiles": "", "generic_scaffold": "", "is_cyclic": False, "num_rings": 0}
    try:
        from rdkit import Chem
        from rdkit.Chem.Scaffolds import MurckoScaffold
        mol = Chem.MolFromSmiles(smiles.strip())
        if mol is None:
            return {"scaffold_smiles": "", "generic_scaffold": "", "is_cyclic": False, "num_rings": 0}
        n_rings = int(mol.GetRingInfo().NumRings())
        is_cyclic = n_rings > 0
        scaff_smi = ""
        gen_smi = ""
        if is_cyclic:
            scaff = MurckoScaffold.GetScaffoldForMol(mol)
            if scaff is not None:
                scaff_smi = Chem.MolToSmiles(scaff)
                try:
                    gen_mol = MurckoScaffold.MakeScaffoldGeneric(scaff)
                    if gen_mol is not None:
                        gen_smi = Chem.MolToSmiles(gen_mol)
                except Exception:
                    gen_smi = ""
        return {
            "scaffold_smiles": scaff_smi,
            "generic_scaffold": gen_smi,
            "is_cyclic": is_cyclic,
            "num_rings": n_rings
        }
    except Exception:
        return {"scaffold_smiles": "", "generic_scaffold": "", "is_cyclic": False, "num_rings": 0}


def compute_scaffold_batch(items: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
    """Batch extraction of Murcko scaffolds from a list of SMILES items."""
    results = []
    for it in items:
        cid = it.get("id", it.get("name", "ITEM_%d" % (len(results) + 1)))
        smi = it.get("smiles", it.get("SMILES", ""))
        info = extract_murcko_scaffold(smi)
        results.append({
            "id": cid,
            "smiles": smi,
            "scaffold_smiles": info["scaffold_smiles"],
            "generic_scaffold": info["generic_scaffold"],
            "is_cyclic": info["is_cyclic"],
            "num_rings": info["num_rings"]
        })
    return results


# --------------------------------------------------------------------------
# 6. Command Line Interface (CLI)
# --------------------------------------------------------------------------
def main():
    parser = argparse.ArgumentParser(description="IDSL_MINT Joint 2215-bit Fingerprint, Embedding & Toxicity Inference Engine")
    parser.add_argument("--mode", type=str, default="predict", choices=["predict", "embed", "smiles_fp", "scaffold"],
                        help="Operating mode: predict (fingerprint + tox), embed (256-d spectrum vector), smiles_fp (candidate structure fingerprint), scaffold (Murcko scaffold extraction)")
    parser.add_argument("--input_json", type=str, required=True, help="Path to input JSON list")
    parser.add_argument("--output_json", type=str, required=True, help="Path to output JSON results")
    parser.add_argument("--weights_path", type=str, default=None, help="Path to PyTorch model weights (.pt)")
    parser.add_argument("--config_path", type=str, default=None, help="Path to model config (.json)")
    parser.add_argument("--tox_model_path", type=str, default=None, help="Optional path to tox_model.joblib")
    parser.add_argument("--include_embedding", action="store_true", help="Include 256-d embedding in predict output")
    parser.add_argument("--device", type=str, default="auto", choices=["auto", "cuda", "cpu"])
    args = parser.parse_args()

    with open(args.input_json, "r", encoding="utf-8") as f:
        inputs = json.load(f)

    if args.mode == "smiles_fp":
        results = compute_smiles_batch(inputs)
    elif args.mode == "scaffold":
        try:
            from rdkit import Chem
            from rdkit.Chem.Scaffolds import MurckoScaffold
        except (ImportError, ModuleNotFoundError) as e:
            sys.exit("Error: RDKit is required for scaffold extraction mode: %s" % e)
        results = compute_scaffold_batch(inputs)
    else:
        if not args.weights_path or not args.config_path:
            sys.exit("Error: --weights_path and --config_path are required for predict and embed modes.")

        pipeline = JointFPInferencePipeline(
            weights_path = args.weights_path,
            config_path = args.config_path,
            tox_model_path = args.tox_model_path,
            device = args.device
        )

        if args.mode == "embed":
            results = pipeline.encode_batch(inputs)
        else:
            results = pipeline.predict_batch(inputs, include_embedding = args.include_embedding)

    with open(args.output_json, "w", encoding="utf-8") as f:
        json.dump(results, f, indent=2)

    print("Successfully processed %d items in %s mode. Output saved to %s" % (len(results), args.mode, args.output_json))


if __name__ == "__main__":
    main()

