# NTSminer 0.1.41

## 大规模多源高分辨质谱融合预训练（550,000 光谱基座大模型）与全类环境污染物泛化升级

> **本版主题**：基于多源海量实验高分辨质谱库（**Spectraverse** 47.7万 Orbitrap 谱图 + **mzCloud** 14.4万环境农药全碰撞能谱图，总计 **549,785 张**训练谱图，覆盖 **>44,000 个独立小分子化学骨架**），在远程 GPU 工作站（NVIDIA RTX 5060 8GB）上完成大规模多源融合预训练（方案 2）。
> **跨越式性能跃升与泛化突破**：
> - **验证集相似度跨越式跃升**：在 71,323 张完全独立未见光谱（严格按 InChIKey 零泄漏分组）上，**Tanimoto 谷本相似度从 0.1379 飙升至 0.7130（相对提升 +417.0%）**，困惑度（Perplexity）从 6.15 降至 2.47，交叉熵损失下探至 0.9041；
> - **广谱环境污染物（农药/除草剂/PFAS）高保真特征感知**：典型全氟烷基酸（PFOA）指纹召回率从 17.9% 提升至 67.9%（Tanimoto 提升 +242.7%），三唑类农药杀菌剂（丙环唑）召回率提升近 4 倍；在 50 个真实水样未知物中精准捕获饱和碳链骨架（100%）、氟原子（98%）、卤代脂肪烃（86%）及磺酸/羧基极性头基；
> - **生产级资产安全与零破坏兼容**：权重资产固化为 FP16 半精度（88.55 MB，严格 <100MB GitHub 限制），CPU 解绑安全加载；`joint_fp_model_config.json` 升级至 1.1.0 并双轨保留基座与多源元数据；全包现有 215 个导出 API 100% 保持完全兼容，核心 228 项单元测试全部满分通过。

---

# NTSminer 0.1.40

## logIE 电离效率预测模型全面优化：特征空间重构、过拟合根治与保留时间（RT）深度融合

> **本版主题**：对内置 **logIE（对数电离效率）预测模型** 与半定量评估流水线进行全面重构与优化。针对此前模型存在的维度灾难（4617 个特征对几百个样本导致的过度过拟合）、R 4.4 下 C++ 指针反序列化崩溃（`XGBAltrepPointerClass`）以及保留时间填 0 缺陷进行系统性根除。
> **保留时间（RT）深度赋能**：在特征重要性分析中，RT 在正模式下 Gain 达 5.5%（居全特征第 3 位），直接表征洗脱时的溶剂化自由能；新增 `gradient_profile` 支持由洗脱程序内插有机相比例。
> **模型瘦身与抗过拟合显著跃升**：
> - **POS 正离子模型**：特征空间由 4617 维大幅缩减至 190 维（166 位 MACCS 指纹 + 21 项物理化学描述符 + RT + MW + Organic_Percentage），Train-Test $R^2$ 差距由 0.360 降低至 0.257（降低 28.6%）；
> - **NEG 负离子模型**：特征空间由 4616 维精简至 355 维（166 位 MACCS 指纹 + 187 项精选 RCDK 描述符 + RT_avg + MW），保留高贡献的含氟特征位（如 `khs.sF`）和拓扑极性参数，Train-Test $R^2$ 差距由 0.283 降低至 0.173（降低 38.9%），测试集 $R^2$ 达 0.610，浓度回算 3 倍误差命中率达 67.5%；
> - **底层架构加固**：通过 `xgboost::xgb.Booster.complete()` 固化模型结构，彻底消除序列化悬空指针；测试套件 1505 个断言 100% 全绿通过。

### 核心改进与模型更新

- **模型特征空间瘦身与重构**：
  - 正离子模式下剔除上千个稀疏指纹位与冗余拓扑指标，精选与 ESI 电离机理紧密相关的物理化学描述符（分子量 `MW`、极性表面积 `TopoPSA`、电荷与折射度 `BCUT` 系列、氢键受体/供体及分子多态性参数）；
  - 负离子模式下保留 PFAS 特征含氟碎片位（`khs.sF` 相对增益达 7.6%）与连接性指数，大幅消除共线性噪音；
  - 重新训练并完整同步更新 `inst/extdata/models/logIE/` 与 `data/logIE_models/` 目录下的正负离子模型及对应 `*_config.json` 规范配置文件。
- **实测保留时间（RT）与梯度深度融合**：
  - `predict_logIE()` 新增 `gradient_profile` 参数，支持通过 `predict_organic_percentage()` 自动根据梯度表内插洗脱有机相比例；
  - `predict_rie()` 自动从输入的特征表或 `NTSProject` 提取实测 `rt` / `RT` / `rt_avg` 与 `organic_pct`，传递给预测引擎；
  - 修复 `predict_rie()` 模型自动加载逻辑，仅在极性明确（`polarity != "auto"` 或从工程中自动识别）时自动加载模型，未指定极性时优雅降级为基于 LogP 的线性近似，100% 保持旧工程与纯数据框单测兼容。
- **底层序列化与架构健壮性加固**：
  - 在 `build_logIE_model()` 与 `save_logIE_model()` 中原生引入 `xgboost::xgb.Booster.complete()`，防止跨版本 R / XGBoost 环境下反序列化丢失外部指针；
  - `build_logIE_model()` 修复极性硬编码为 `"positive"` 的历史缺陷，支持 `polarity = c("positive", "negative")`，并新增 `metadata_matrix` 支持传递外部色谱元数据；
  - 保留自定义调优参数 `tune_params$objective`（如 Huber 鲁棒损失函数 `reg:pseudohubererror`），不再被硬性覆盖为 MSE。

---

# NTSminer 0.1.39

## 正向高分辨 MS/MS 质谱图全流程预测引擎（FraG-FIORA 范式）与 L2b 级计算机辅助注释

> **本版主题**：新增**正向高分辨 MS/MS 质谱图全流程预测引擎（Forward MS/MS Spectrum Prediction，基于 FraG-FIORA 混合范式）**。针对非靶向筛选中大量缺乏参考标准品实测谱图的候选污染物，直接从化合物化学结构（SMILES）正向推导其在不同碰撞能量（NCE 10, 20, 40 或 Merged 合成谱）下的完整高分辨 MS/MS 质谱图（精确到小数点后 4 位的单同位素 $m/z$、相对丰度及断裂机理/中性丢失标签），支撑候选物正向打分重排与 Schymanski Level 2b（计算机模拟/预测库匹配）级注释。
> **广谱污染物通用性与 PFAS 深度微调**：
> - **广谱污染物支持**：全面覆盖各类环境新污染物（喹诺酮类、磺胺类、四环素类、农药、有机磷酸酯 OPEs、邻苯二甲酸酯 PAEs 及其他复杂芳香/杂环小分子）；
> - **实测 In-House PFAS 库物理微调**：基于 73 种标准品实测高分辨谱库（`NTS_InHouse_PFAS.rds`），吸取此前 CFM-ID 微调经验，采用物理约束 L-BFGS-B 优化器针对含氟特征断裂（C-F 链阶梯断裂、-CO2、-SO3、-HF 消除、母离子存活率）进行微调训练，严格维持多碎片信号召回率 $>50\%$，并输出专用权重固化于 `inst/extdata/models/forward_ms2/`。
> **向后兼容与平滑扩展**：全包原有 211 个导出函数 100% 保持兼容，导出函数集由 211 增至 215（新增 `predict_ms2_spectrum`、`rescore_candidates_by_predicted_spectrum`、`annotate_by_insilico_prediction`、`plot_predicted_mirror`）；S4 `NTSProject` 槽位规范与数据流严格兼容。

### 核心功能与架构升级

- **`predict_ms2_spectrum()` 核心正向预测函数**：
  - **图论组合剪枝 + 深度倾向性赋权**：吸收 FraGNNet 严格质量守恒与组合子图分割优势，根除经验模型常见的非物理假峰（Zero Hallucination）；同时吸收 FIORA 局部化学微环境断裂倾向性与气相电荷竞争模型；
  - **自动化学路由（`model = "auto"`）**：针对氟原子数 $\ge 3$ 的全氟/多氟化合物自动启用 in-house 微调模型权重（`pfas_inhouse`），其他通用小分子启用通用权重（`general`）；
  - **能量梯度调控**：支持输出 NCE 10、NCE 20、NCE 40 单一能量谱，或加权融合的 `"Merged"` 全能量合成谱；
  - **面向对象的专用 S3 输出**：单化合物输出 `nts_predicted_spectrum`（继承 `data.frame`），多化合物输出具名列表，配备格式化终端打印方法 `print.nts_predicted_spectrum()`；
  - **优雅降级保障**：在无 Python 或未安装 RDKit 环境下发出友好警示并返回结构完备的空谱图对象，避免硬崩溃。
- **`rescore_candidates_by_predicted_spectrum()` 候选物谱图匹配与重排序**：
  - 将实测 MS2 碎片谱图与候选化合物的正向预测谱图进行高通量贪婪匹配（容差 `mz_tol = 0.02` Da）；
  - 输出平方根加权余弦相似度（`cosine_sqrt`，降低动态范围主导、强化微弱特征碎片匹配）、线性余弦（`cosine_linear`）、解释强度占比（`explained_intensity`）与综合打分（`score_pred`），输出规范化降序排名长表。
- **`annotate_by_insilico_prediction()` 自动化工程级 L2b 注释流水线**：
  - **贯通 MS1 与 MS2 证据链**：对 `NTSProject` 内未知特征首先通过前体离子质量在 10 ppm 容差下初筛候选结构清单（默认接入包内 2855 条三库去重 PFAS suspect 清单）；
  - **动态预测与正向过滤**：对拥有实测 MS2 的特征现场调用 FraG-FIORA 生成各候选物的高分辨 MS/MS 谱图并快速计算匹配指标，按阈值（`min_cosine >= 0.30` 或 `min_explained_int >= 0.25`）保留 Top-N 候选；
  - **Schymanski L2b 级注释写入**：以 `db_type = "prediction"` 优雅回写至 `@db_annotations`，天然接入包内置信度评估体系（`assign_confidence_levels()` 自动评定为 L2b），并详细记入工程日志 `@log`；
  - **既有注释保护**：`skip_annotated = TRUE` 严格保护已获得 L1（in-house）或 L2a（online）确证的特征，杜绝低置信度预测覆盖高置信度实验命中。
- **`plot_predicted_mirror()` 镜像对比图可视化**：
  - 生成高质量 `ggplot2` 镜像对比图，上方实测谱（红），下方 FraG-FIORA 预测谱（蓝），并在峰尖旁竖排标注匹配的单同位素 $m/z$，副标题实时呈现 Cosine (sqrt)、解释强度及匹配峰统计。
- **Python 推理引擎扩展（`inst/python/frag_fiora_engine.py`）**：
  - 原生支持单 SMILES 模式与 `--batch_json` 高吞吐批处理模式，避免频繁创建 Python 进程的 Windows 进程开销（批处理速度提升数十倍）；内置参数缓存机制，支持单批混合不同模型权重。

---

# NTSminer 0.1.38

## 无库未知物类似物与化学骨架自动推导（DreaMS 范式）：语义空间邻居检索与骨架修饰解析

> **本版主题**：新增**无库未知物类似物与化学骨架自动推导（Analog & Scaffold Retrieval）**模块。针对非靶向质谱中绝大多数无谱库直接命中、处于 Schymanski Level 4/5 级的未知新污染物，借助 Transformer 编码器抽离的 256 维高维谱图语义连续表征空间，检索其语义最近邻确证类似物，结合精确质量差（$\Delta m/z$）与通用环境化学转化字典，提取母核骨架（Bemis-Murcko Scaffold）及化学分类标签，打破未知物无信息孤立瓶颈。
> **通用小分子分析定位**：适用于各类环境小分子污染物（PFAS、农药、医药品与个人护理品 PPCPs、工业助剂、增塑剂及其环境转化产物），内置 28 种通用环境相 I/II 代谢与同系物化学修饰转化规则。
> **向后兼容与平滑扩展**：全包原有 208 个导出函数 100% 保持兼容，导出函数集由 208 增至 211（新增 `default_transformation_rules`、`infer_scaffold_from_semantic_neighbors`、`get_inferred_scaffolds`）；S4 `NTSProject` 槽位规范严格保持不变。

### 核心功能与架构升级

- **`infer_scaffold_from_semantic_neighbors()` 核心推导函数**：
  - **输入多态性**：支持 `NTSProject`、特征数据框或列表等多态输入；支持直接复用工程内预计算嵌入或外部自定义嵌入矩阵；
  - **多源参考类似物对齐**：支持自动识别工程内已确证参考种子（`@db_annotations` 中携带有效 SMILES 的化合物）、用户指定种子集（`anchor_features` / `anchors` 表），以及外部参考数据库（`reference_db` 支持 `"inhouse"`、`"online"` 或自定义 `NTSDatabase`）；
  - **极速稠密检索与 Top-K 过滤**：利用底层 BLAS 矩阵乘法，在毫秒级内完成数千未知特征与参考候选集的余弦距离比对，严格按 `min_semantic_score`（默认 0.70）与 `top_k`（默认 3L）截断；
  - **通用环境转化与精确质量差匹配**：内置 28 项环境转化规则（包含羟基化、脱羧、脱氟、硫酸化、葡萄糖醛酸化、±CF2/CH2 同系物序列步进、醚键断裂等），精确计算 $\Delta m/z$ 与 ppm 误差，自动判定转化修饰类型；
  - **化学母核与分类推导**：Python 推理引擎扩展 `--mode scaffold`，通过 RDKit 提取 Bemis-Murcko Scaffold（环系骨架）与通用拓扑（Generic framework），在无 Python 环境下优雅降级为纯 R 环系特征识别，结合化学分类提供结构解析依据；
  - **结构化结果合成与工程回写**：生成人类可读的综合推导摘要（例如 `"Analog to PFOA (sim=0.980); shift: -CF2 (Δmz=-49.9968, 0.4 ppm); Acyclic; Class: Perfluoroalkyl carboxylic acids"`）；`write_to_project = TRUE` 时将主推导表与详细候选表优雅写入 `@physchem_proxy` 并记入工程日志 `@log`。
- **便捷查询器 `get_inferred_scaffolds()`**：
  - 一键从 `NTSProject` 中提取推导主表（`type = "summary"`）或全量 Top-K 候选匹配长表（`type = "analogs"`）。
- **Python 推理引擎升级（`inst/python/joint_fp_engine.py`）**：
  - 新增 `extract_murcko_scaffold` 与 `compute_scaffold_batch`，CLI 支持 `--mode scaffold`，实现高通量批处理结构母核提取，缺失 RDKit 时输出标准错误保障 R 侧优雅降级。

### 独立审计与防御加固

- **经资深独立子 Agent 多维度穿透审查**，通过严格把关验收并落地多项防御加固：
  - **[P1] 修复降级锚点特征排序错位**：使用 `match()` 严格对齐 `af_ids` 与特征表行序，彻底规避因顺序错乱导致的化合物身份误配；
  - **[P1] 防御 `delta_mz` 非标量崩溃**：增加标量与有限性严格校验，杜绝 `numeric(0)` 触发条件判定异常；
  - **[P2] 未命名嵌入矩阵行名自动推导**：支持无行名嵌入矩阵在长度对齐时自动绑定特征 ID，防止锚点被意外清空；
  - **[P2] 非有限数值（`NA/NaN/Inf`）防御性剔除**：在矩阵归一化前过滤非法数值，防止 BLAS 乘法污染；
  - **[P2] 严格锚点去重与自匹配排除**：避免已知特征在对照测试中与自身产生同名无意义匹配；
  - **[P2] Python RDKit 降级保障**：在 `--mode scaffold` 中增加显式模块检测，确保在缺少 RDKit 的环境中平滑回退纯 R 正则检测；
  - **测试保障**：新增 `tests/testthat/test-semantic-scaffold-inference.R`（72 个断言 100% 通过），全包 5 套语义网络核心单测（228 个断言）全绿通过。

### 真实非靶向环境 PFAS 质谱工程评测

- 基于真实非靶向高分辨项目 `project_final.rds`（4893 MS1 特征、966 纯化 MS2、114 库确证种子）进行全量实测：
  - **未知物类似物推导覆盖率 100%**：852 个无库命中的未知 MS2 特征在 2.05 秒内**全部成功匹配到高相似度确证类似物（相似度 $\ge 0.70$）**，输出 2,556 条 Top-3 候选关系；
  - **母核拓扑分布**：770 个特征识别为无环链状骨架（全氟烷基链衍生），82 个特征推导出含环系/芳香母核；
  - **标杆案例确证**：83 个未知物识别为 PFOA 类似物（余弦相似度达 0.9910），94 个未知物识别为 PFOS / EtFOSA 类似物（余弦相似度达 0.9899），成功捕获同系物 $\pm\text{CF}_2$、脱羧 $-\text{CO}_2$ 与同分异构体。

---

# NTSminer 0.1.37

## 混合分子网络（Hybrid Molecular Network）：硬规则边与深度质谱语义嵌入边无缝融合

> **本版主题**：新增**混合分子网络（Hybrid Molecular Network）**架构。将基于硬规则的经典分子网络（共有碎片 DFI、中性丢失 NL、谱内质量差 MD、同系物序列 KMD、生物/非生物转化关系 Transformation、源内碎片 ISF）与基于 Transformer / 深度表征模型的质谱语义嵌入边（256 维稠密特征向量余弦相似度）无缝融合。
> **通用小分子分析定位**：适用于 PFAS、农药、医药品与个人护理品 (PPCPs)、增塑剂、工业添加剂等各类环境小分子污染物，算法与接口严禁硬编码限定在单一类别。
> **向后兼容与优雅降级**：原有 208 个导出函数签名与位置参数 100% 保持兼容；`include_semantic = FALSE` 默认关闭；原有网络分析脚本行为逐字节不变。

### 核心功能与架构升级

- **`get_network_data()` 与 `generate_network_data()` 语义层打通**：
  - 新增参数 `include_semantic = FALSE`、`min_semantic_score = 0.70`、`semantic_top_k = 3L`、`embeddings = NULL`、`python_path = NULL`；
  - `type` 切片新增 `"semantic"` 选项，支持单独提取语义关联子网络；
  - 支持直接复用工程 `@physchem_proxy$ms2_embeddings` 或传入外部预计算嵌入矩阵（如 DreaMS 1024 维或 Spec2Vec 300 维），利用底层 BLAS 矩阵乘法极速完成全特征相似度矩阵计算（$O(N^2 D)$）；
  - 遵循 GNPS Top-$K$ 修剪准则，每个节点只保留最强的前 $K$ 个语义邻居，防止局部高相似度簇导致图边爆炸；
  - 当项目缺少 `@clean_ms2` 时，语义边构建自动解耦、独立成边，并在混合模式下优雅降级（发出警告并跳过峰边，语义边照常构建）。
- **多层网络边折叠与证据融合（`.nts_merged_network_edge_collapse()`）**：
  - 当节点对同时具有硬规则物理边（如共有碎片或同系物）和 MS2 语义边时，自动折叠为一条复合边，保留各自的证据列（`shared_fragments_*`、`shared_nls_*`、`semantic_similarity`），边打分 `edge_score` 自动取各层最大值，实现多重证据交叉强化。
- **可视化层深度集成（`plot_nts_network()` / `plot_molecular_network()`）**：
  - 语义边赋予专属的中紫配色（`#9370DB`），线型设为辨识度高的虚线（`c(5, 5)`），箭头呈现无向线段（`arrows = ""`）；
  - 边悬停提示框（Tooltip）动态渲染 3 位小数的语义余弦相似度（`Semantic similarity: 0.xxx`）；
  - `plot_nts_network()` 完整透传语义参数，无缝衔接现有交互式图谱与子网络提取功能（`extract_subnetwork()`）。

### 独立审查与测试保障

- **通过资深独立子 Agent 多维度穿透审查**：审查覆盖向后兼容性、领域普适性、数学计算严谨性、边界防御与性能内存；
- **边界防御加固**：
  - 修复 partial `@ms2_annotations` 时语义边端点被过度收窄的潜在缺陷；
  - 增加未命名嵌入矩阵的自动行名匹配与非有限数值（`NA/NaN/Inf`）防御性剔除；
  - 防御非法或 `NA` 的相似度阈值输入。
- **单测覆盖**：
  - 新增 `tests/testthat/test-hybrid-molecular-network.R`（53 个断言 100% 通过）；
  - 既有网络测试 `test-network-s1-s3.R`（53 个断言）全绿通过。

---

# NTSminer 0.1.36

## 候选结构指纹重排序与质谱语义嵌入（Spec2Vec 范式）

> **本版主题**：新增 **MS2 联合指纹驱动候选库反向重排（Candidate Re-ranking）** 与 **基于 Transformer 稠密特征的质谱语义相似度嵌入（Spec2Vec 范式）**。
> 无任何破坏性变更；原有 204 个导出函数全量兼容；**导出函数集由 204 增至 208**。

### 新增导出 API

- **`rescore_candidates_by_joint_fp(candidates, predicted_fp = NULL, smiles_col = "SMILES", weight_fp = 0.30, score_col = "score_total", ...)`**：
  将 Transformer 预测的 MS2 联合分子指纹（2215 位）与候选库化学结构（SMILES）进行位级谷本相似度（Tanimoto Similarity）匹配，并将该相似度加权融合原有谱库得分或打分（`score_total`），生成综合再打分 `score_reranked` 并完成特征内候选重新排序。
  - 支持 `NTSProject`（直接从 `@physchem_proxy` 读取预测指纹）、数据框或具名位向量列表；
  - 候选指纹计算自动降级：候选数据框直接提供位列 -> Python RDKit 批量指纹 -> `rcdk` 本地回退；
  - 针对高相似度（$\ge 0.70$）候选自动标记 `FP_MATCH`。
- **`compute_candidate_joint_fp(smiles, python_path = NULL, verbose = FALSE)`**：
  为候选化学结构批量计算标准的 2215 位联合二进制指纹（0..2047 位 ECFP2048 + 2048..2214 位 MACCS167），输出激活位索引列表。
- **`extract_ms2_embeddings(object, precursor_mz = NULL, peaks = NULL, ...)`**：
  基于 Transformer 编码器抽离未知物质谱图的 256 维稠密特征向量（Global Mean Pooling + L2 Normalization）。单谱仅需单次前向编码（无需经过自回归解码器），推理极快（单谱 < 2 ms）。支持多态输入并可直接写回 `NTSProject` 的 `@physchem_proxy` 槽位。
- **`match_spectral_semantic(query, reference, top_k = 5L, min_similarity = 0.0, ...)`**：
  利用底层 BLAS 矩阵乘法（内积等价于 L2 归一化后的余弦相似度），在数十毫秒内实现数千未知物与数千参考谱的语义相似度检索，突破传统硬点积峰对齐算法对仪器质量偏差和微小碎裂差异的敏感限制。

### 注释排序与平票仲裁引擎打通

- **`rank_annotation_candidates()`**：
  新增 `predicted_fp = NULL` 与 `weight_fp = 0.25` 参数。在多源候选合并（实测库、预测库、MetFrag、网络推断）中，根据指纹相似度自动打破同分异构体及相同置信等级候选间的平票（Tie-breaking），在 `candidates` 与 `final` 决策表中输出 `fp_similarity`、`score_reranked` 与 `fp_sim_best`。
- **`merge_annotation_sources()`**：
  透传 `predicted_fp` 与 `weight_fp`，一键式完成四路证据融合与深度学习指纹重排序。

### Python 推理引擎扩展（`inst/python/joint_fp_engine.py`）

- 新增 `MS2FP_Model.encode_spectrum()` 与流水线 `encode_spectrum()` / `encode_batch()`；
- 新增 `compute_smiles_joint_fp()` 与批量转换支持；
- CLI 支持 `--mode {predict, embed, smiles_fp}` 与 `--include_embedding` 选项，100% 保持旧命令参数兼容。

### 单元测试与验证

- 新增 `tests/testthat/test-joint-fp-reranking.R`（25 个断言 100% 通过）；
- 新增 `tests/testthat/test-spectral-semantic.R`（28 个断言 100% 通过）；
- 回归验证 `test-joint-fp-prediction.R`（50 个断言）与 `test-annotation-rank.R`（102 个断言）全绿通过。

---

# NTSminer 0.1.35

## 深度学习 MS2 联合分子指纹与毒性预测模块（IDSL_MINT Joint Model）

> **本版主题**：新增**深度学习 MS2 质谱图直推 2215 位联合分子指纹（ECFP2048 + MACCS167 偏移 2048 位）及下游 QSAR 环境毒性风险初筛引擎**。
> 无任何破坏性变更；原有 201 个导出函数及其打分算法逐值不变；**导出函数集由 201 增至 204**。

### 新增导出 API

- **`load_joint_fp_model_config(config_path = NULL)`**：
  纯 R 读取内置或用户指定的 `joint_fp_model_config.json`（**无需 Python / PyTorch**），返回 S3 对象 `nts_joint_fp_config`。
  打印展示包含网络超参数（148,003 质荷比分箱、2,218 词表、256 隐藏维、8头4层 Transformer）、5 轮全量训练收敛曲线（Loss 7.9361 -> 1.3329，留出集 Tanimoto 0.6802），以及在 **174 个完全隔离的 PFAS 未见骨架（3,858 条实测谱图）** 黄金基准集上的独立评测指标与下游毒性拟合 Spearman 秩相关系数（$\rho = 0.7303$）。
- **`predict_joint_fingerprint(object, precursor_mz = NULL, peaks = NULL, ...)`**：
  将未知物原始 MS2 碎片谱图与母离子 $m/z$ 经自包含 Transformer 贪婪自回归解码为 2215 位联合指纹（ECFP 0..2047 位，MACCS 2048..2214 位），自动解耦出关键环境官能团（氟原子、脂肪族 C-F 键、全氟烷基链、羧酸根、磺酸根等），并联动下游 Random Forest QSAR 模型预测毒性分值（0.0~1.0）与风险档次（High / Moderate / Low）。
  - **多态输入支持**：支持 `NTSProject`（自动提取 `@clean_ms2` 或 `@ms2data`，结果可回写 `@physchem_proxy` 并更新 `@log`）、二维峰表矩阵/数据框、多谱列表。
  - **优雅降级机制**：自动探测 Python / PyTorch（优先级：入参 `python_path` > `NTSMINER_PYTHON` > `RETICULATE_PYTHON` > 系统 PATH），环境缺失时给出明确指引，严禁硬崩溃。
- **`predict_toxicity_from_ms2(object, ...)`**：
  面向环境非靶向未知物高通量初筛的高层便捷封装，一步完成从 MS2 质谱图到毒性风险分级与官能团识别。

### 预训练深度学习资产固化（Option A）

- **随包内置 3 大核心资产（`inst/extdata/models/joint_fp/`）**：
  - `joint_fp_model_config.json`：网络超参数、词表映射与黄金基准评测元数据。
  - `joint_fp_model_weights_fp16.pt`：纯推理半精度权重（88.55 MB，小于 GitHub 100MB 限制）。
  - `tox_model.joblib`：下游环境毒性 QSAR 分类模型（6.23 MB）。
- **自包含轻量 Python 推理引擎（`inst/python/joint_fp_engine.py`）**：
  - 纯 PyTorch 实现，零外部依赖，集成因果掩码与 `expandable_segments` 显存防护。
  - 支持单谱与批处理 CLI（`--input_json` / `--output_json`），进程隔离，调用后立即释放 GPU/系统内存。

### 单元测试与兼容性

- 新增 `tests/testthat/test-joint-fp-prediction.R`（50 个断言 100% 通过）。
- 全包原有 204 个导出函数全量回归（1,263 个测试通过），完全保持向后兼容。

---

# NTSminer 0.1.34

## SIRIUS 模块修复（元素参数语义 / CANOPUS 列名 / 外部进程失败被吞 / output_dir 误删）

> **本版主题**：修 `docs/analysis/NTSminer_SIRIUS模块审阅_0.1.33.md` 的阶段 A 11 项。
> **有行为变化**（见「破坏性 / 行为变化」）；`@sirius_fingerprints` 的 14 列输出契约
> **不变**，不新增、不删除任何列。

### 破坏性 / 行为变化

- **`elements` 默认值与语义（A-P0-1）**：`run_sirius()` / `compute_sirius_fingerprints()`
  的 `elements` 默认值由 `"CHNOPFClBr"` 改为 **`"SClBrSeI"`**。语义明确为
  **`--elements-considered` 的"稀有元素候选集"**——SIRIUS 始终考虑 C/H/N/O/P，
  该参数只负责**再放开哪些可由软件自动判定的稀有元素**。旧默认**不含 S**，
  PFAS 磺酸/磺酰胺类会被漏判；新默认补上 S 并追加 Se / I，C/H/N/O/P 不再显式列出。
  `NULL` → 直接**省略该参数**（用 SIRIUS 默认）；传入的符号串里若**不含 S** 会
  **告警**（提醒 PFAS 场景下 S 通常是必需的）。元素符号**区分大小写**（`Cl`/`Br`/`Se` 不做大写化）。
  > ⚠️ 该参数在 SIRIUS CLI 中的**确切语义**尚未在本机实测（沙箱跑不了 Windows
  > `sirius.exe`）。本版按"保守实现乙"落地，**实测待补**，详见修复证据文件。
- **外部进程失败不再被吞（S-P1-1 / S-P1-4）**：此前 `run_sirius()` 的 `system2()`
  返回值**直接丢弃**，SIRIUS 崩了也当成功。现在：
  - `run_sirius()`：非零退出码 **`stop()`** 并给出命令行与日志路径；
  - `compute_sirius_fingerprints()`：非零退出码 **告警并记录到 `batch_status`**，
    不中断整批（可继续解析其它批次的产物）；
  - 两者都把 SIRIUS 的 stdout/stderr **落盘**到 `output_dir`（`sirius_cli.log` /
    `sirius_batch_<i>.log`），不再在脚本化运行时丢弃（S-P1-2）。
- **`output_dir` 不再无条件删除（S-P1-3）**：此前 `unlink(output_dir, recursive = TRUE)`
  对任何路径都执行。现在只清理由 NTSminer 自己创建的目录（以 `.nts_sirius_output`
  标记文件识别）；**已存在但不是 NTSminer 创建的目录会报错**，除非显式传
  **`overwrite = TRUE`**（此时告警后仍清理）。`run_sirius()` 与
  `compute_sirius_fingerprints()` 均新增 `overwrite = FALSE` 入参。

### 结果变化（数值填充口径）

- **`explained_peaks` 改为整数峰数（A-P1-3）**：此前取的是 `explainedIntensity`
  （0–1 的**比例**，还被 `as.integer()` 截成 0/1），现在取 **`numExplainedPeaks`**
  （SIRIUS 6 `formula_identifications.tsv` 第 10 列，**整数**峰数）。该列的值会从
  0/1 变成真实计数。
- **CANOPUS 四列在 SIRIUS 6 下可填（A-P1-2）**：NPC 列解析改为按真实表头
  （`NPC#pathway` / `NPC#superclass` / `NPC#class`，均带 `Probability` 并列列）**精确匹配**，
  拼接格式为 `pathway; superclass; class`。修复了此前**双重转义**的正则导致的列名匹配失败。
  ⚠️ SIRIUS 6 的 `canopus_formula_summary.tsv` **没有 `NPC#subclass` 列**，故
  `canopus_subclass` 在这两版 SIRIUS 下**始终为 `NA`**——这是**预期行为**，不是缺陷
  （`ClassyFire#subclass` 是另一套分类体系，**不会**被拿来填 `canopus_subclass`）。
- **`confidence` 显式告警（A-P2-8）**：SIRIUS 6 的 **formula** 输出**无任何 confidence 列**，
  `confidence` 将为 `NA`；现在会在解析时**明确告警**（而不是静默留 NA）。若某版本
  确实提供了 `ConfidenceScoreExact` / `ConfidenceScoreApproximate`，优先取**后者**。

### 修复（不改变既有调用）

- **A-P1-5 API 指纹按候选行回填**：`extract_sirius_fingerprints_api()` 回填指纹时按
  返回候选的 `rank` / `molecularFormula` / `adduct` **定位到确切行**，不再固定写
  `[1L]`（此前会把指纹写错行；命中多行时告警"仅填首行"）。
- **A-P1-4 CANOPUS 合并按特征键**：合并键加入 `mappingFeatureId`（经 `.nts_sirius_repair_id()`
  规范化），**同分子式不同特征不再串味**；表内无 `mappingFeatureId` 时回退并告警。
- **S-P1-5 `batch_size` 校验**：`0` / 负数 / `NA` / `Inf` / 字符串 / 长度 > 1 一律**报错**；
  `NULL` 与非整数正数（如 `20.0`）正常接受。
- **`rbindlist(fill = TRUE)`（A-P2 项）**：批次元数据由 `do.call(rbind, ...)` 改为
  `rbindlist(fill = TRUE)`，各批列集不一致时不再报错，缺失列按要求补齐为 14 列。
- 新增内部助手便于测试与复用：`.nts_sirius_run_cli()` / `.nts_system2_capture()` /
  `.nts_sirius_prepare_output_dir()` / `.nts_sirius_read_canopus()` / `.nts_sirius_col_map()`。
  新增 **`tests/testthat/test-sirius-parsing.R`**（夹具为真实 SIRIUS 6.3.4 输出，
  **不依赖真实 SIRIUS / 不起服务 / 不联网**）。

# NTSminer 0.1.33

## 包内 xgboost 模型改为新版序列化格式（此前现代 xgboost 无法加载）

> 包内 `inst/extdata/models/` 下的 8 个 xgboost 模型此前存的是 **xgboost 1.7.9 旧式
> UBJSON**（顶层 `{"Config","Model"}`）。这种布局 **xgboost ≥ 2.0 只认统一布局
> `{"learner":{...}}`**，加载会报 `'xgb.Booster' object is corrupted or is from an
> incompatible XGBoost version`；**xgboost 1.7.6 同样不认**（PyPI 也没有 1.7.9）。
> 也就是说：**在一个装了现代 xgboost 的环境里，`predict_rt()` / `predict_logIE()`
> 此前根本跑不通**——不是算法问题，是模型文件本身读不进来。
>
> **新格式由 xgboost 3.4.2 写出，需要 xgboost ≥ 2.0。**

- **修复**：`inst/extdata/models/{rt,logIE}/` 下 8 个模型（neg / pos / xgb / `_rcdk` /
  `_old` 变体）全部重存为统一布局。做法：读旧式 UBJSON → 取 `Model.learner` 并用
  `Config.learner` 补缺键 → 输出统一 JSON → `xgb.load()` → `saveRDS()`。
  `load_rt_model()` / `load_logIE_model()` 的返回结构与 `model_version` 均不变。
- **数值不变**：用**独立的 Python 树遍历实现**对重存前（旧 UBJSON 原始字节）与重存后
  （新 rds，xgboost 预测）逐模型比对 3 组输入，**最大相对差 1.6e-6**（float32 量级）。
  另用真实 est 数据 A/B：同一份结果表**逐单元格完全一致**（324 候选 / 137 特征）。
- **副作用**：模型文件在现代 xgboost 下可直接加载，不再需要调用方自建转换层。
  若环境仍是 xgboost 1.7.x，则需要升级到 ≥ 2.0 才能加载（升级前无论如何都读不了）。


## RT / logIE 预测模块审阅修复（残留泄漏 / 默认模型缺列 / 元数据缺失 / 口径不一致）

> **本版主题**：修 `docs/analysis/NTSminer_RT与logIE预测审阅_0.1.32.md` 列出的问题。
> **有行为变化**：模型评估指标与 RIE 预测值会变（见「破坏性 / 行为变化」）；其余为
> 新增字段与诊断信息，**不删列、不改既有参数默认**。

### 破坏性 / 行为变化（数值会变）

- **P1-1（残留泄漏）**：`build_rt_model()` / `build_logIE_model()` 改为
  **先切分、后拟合特征变换**——近零方差 + 高相关筛选与中心化·标准化**只在训练行上
  学习**，其余行仅被*套用*。此前变换在切分前拟合，测试行的分布参与了缩放常数估计，
  测试指标偏乐观；修复后 **test R²/RMSE 以及 `predict_rt()` / `predict_logIE()` 的
  预测值都会变化**。新增内部助手 `.nts_rt_apply_preproc(mat, prep)` 供复用。
- **P0-1（默认 logIE 模型缺列）**：包内 logIE 模型需要 `RT`/`RT_avg`、`MW`、
  `Organic_Percentage` 等元数据列，而唯一包内调用路径此前只传 SMILES，缺列被静默
  **填 0**。现在 `predict_logIE()` **新增 `rt` / `mw` / `organic_pct` 入参**
  （`RT`/`RT_avg` 取 `rt`，`MW` 取 `mw`）；`MW` 未给时由结构**推导**（`rcdk` 描述符，
  纯结构量，**绝不填 0**）；仍缺的列会**明确列出**并说明估算方式。同时 logIE 路径的
  `.nts_compute_descriptors()` 一律传 `prune = FALSE`——该筛选按**批**评估，会在少
  样本批里把训练时用到的列删掉再被静默填充（**这才是"缺列"的真正根因**）。
  **`predict_rt()` 的缺列填充**同样由原始 0 改为**训练集中心**（标准化后即 z=0 的
  中性值；填 0 会变成很大的 z 值）。
- **P2-5（区间外推）**：`predict_rt(return_interval = TRUE)` **只对域内化合物**给
  区间（模型 test RMSE 的 ±2σ），域外行 `lower`/`upper` 置 `NA`（此前一律给区间 =
  虚假精度）。**不带适用域的旧模型仍照旧给区间**（向后兼容）。
- **P2-6**：`predict_logIE()` 的极性不符由**警告**改为**报错**（一个模型只服务一个
  极性，矛盾时不能"照猜不误"，`test-logIE-prediction.R` 同步更新）。

### 新增 / 修复（不改变既有调用）

- **P1-3**：`build_logIE_model()` 新增**特征筛选**（`feature_screen = TRUE` 默认，
  训练行上学习）与**分组切分**（新增 `group` / `cluster_similarity`，与
  `build_rt_model()` 对齐；分组不 informative 时回退按 logIE 分层切分）。此前
  4617 特征 × ~480 样本**无筛选**且用分层（非分组）切分，测试误差偏乐观。模型对象
  新增可追溯字段 `split_strategy` / `n_groups` / `preproc_center` / `preproc_scale`，
  `predict_logIE()` 据此对输入做与训练**一致**的标准化。
- **P1-2（RT 模型元数据）**：`build_rt_model()` 现把 `model_version` 写为
  **`"0.2.0"`**（新增常量 `.NTS_RT_MODEL_VERSION` / `.NTS_RT_MODEL_LEGACY_VERSIONS`），
  并记录 `app_domain` / `split_strategy` / `n_groups`；`print.nts_rt_model()` 新增
  **`App. domain`** 一行（无适用域时提示 "predictions are extrapolation"）；
  `load_rt_model()` 的版本校验改为 **NULL 安全**（缺失版本只警告、不报错）。
- **P2-3**：`evaluate_rt_model()` / `evaluate_logIE_model()` 在预测数与观测数不一致时
  **不再按位置硬拷**（会静默错配），改为置 `NA` 并警告，受影响行排除出指标计算。
- **P2-4**：`annotate_rt_prediction()` 的结果按 **`feature → |rt_delta|`** 排序
  （此前把整个 `@db_annotations` 按 `abs(rt_delta)` **全局**排序，跨特征乱序）。
- **P3-1**：新增内部助手 `.nts_resolve_rt_model_path(pol_abbr)`（`R/08_risk_assess.R`），
  与 logIE 版镜像，解析 `inst/extdata/models/rt/rt_model_inhouse_<pol>.rds`。

### 数据 / 元数据 / 测试

- 包内 **logIE 四个 `*_config.json`** 重写为**核对无误**的内容（特征数、训练/测试样本
  数、test R²/RMSE/MAE、cv R²），并注明 `_old` / `_rcdk` 与规范文件**逐字节相同**
  （md5 `d95e750720aedaaad7b67c2a159df1c2`）；**RT 两个 `*_config.json`** 补 `file`
  字段并重排；新增 **`rt_model_inhouse_xgb_config.json`**（通用 56148 样本模型，注明
  **未记录适用域**，故 `predict_rt()` 对其返回 `NA`）。
- 测试新增 **`tests/testthat/test-rt-logIE-fixes-0133.R`（78 断言）**；
  `test-isotope-scoring.R` 的「空 feature_table」用例补 `skip_if_not_installed("Rdisop")`
  （与其同级 10 个用例一致；此前在无 Rdisop 环境下会先撞上后端报错）。

### 已知残留（本版未改）

- 包内 logIE 模型**未回填** `preproc_center` / `feature_centers`（需重训生成）；带中心
  缺失的模型按**不缩放**预测，与 0.1.32 逐位一致。
- `nAtomLAC`（`rcdk` 3.8.2 未提供）无法生成，仍按 1 列填充并**告警**。
- 由**预计算指纹矩阵**构建的模型不记录 `fp_types`，故 `predict_rt(smiles = ...)` 无法
  为其现算指纹（须传 `fingerprint_matrix`）；`evaluate_rt_model()` 不接受矩阵入参，
  对这类模型不可用——**未改**（超出本次审阅范围）。

### 验证

- 全量 testthat：**1072 通过 / 0 失败 / 0 报错 / 0 警告 / 37 跳过**（42 个文件）。
- 沙箱内以桩 `ranger` 打通 `rf` 后端，实跑核对：切分 → 训练行拟合变换 → 域内/域外
  预测 → 评估 → 存读 round-trip，逐条断言通过。


# NTSminer 0.1.32

## 注释合并 / 优先级排序：多路证据 → 每特征一个最优注释（新增导出 API）

> **本版主题**：把一条龙跑出来的**四路注释证据**（实测谱库命中 / 预测谱库命中 /
> MetFrag 预测 / 网络·同系物传播）合并成 **一张候选排序表 + 每特征一个最优注释**，
> 并给出**人工复核清单**。**只读派生**：不写工程、不写 `ckpt_*.rds`、不写
> `@network_results`。

### 新增

- **`rank_annotation_candidates(candidates, metfrag_hits = NULL, network_layers = NULL,
  feature_meta = NULL, suspect = NULL, verbose = TRUE)`** —— 纯数据框核心，无副作用。
  `candidates` 为库命中候选长表（特征列 `feature` / `name` / `feature_id` + `ref_id`），
  可选 `metfrag_hits` / `network_layers` / `feature_meta` / `suspect`；
  返回四个 `data.frame`：`candidates`（逐特征×候选，含 `rank_in_feature` / `is_top` /
  `base_level` / `decision_tier` / `decision_label` / `n_families` /
  `family_F1` / `family_F2` / `family_F4` / `flags`）、`final`（每特征一行，含
  `tied_alternatives` / `review_reason` / `conflict_status` 与各族首选 `fam1_best` /
  `fam2_best` / `fam4_best` / `lib_best_ref` / `lib_best_scope` / `mf_top1_ref` /
  `mf_vs_lib` / `net_ident_ref`）、`conflicts`（`review_reason` 非空的子集）、
  `summary`（按 `decision_label` 计数）。
- **`merge_annotation_sources(run_dir, out_dir = file.path(run_dir, "merged"),
  suspect_csv = NULL, multilib_set = NULL, write_files = TRUE, verbose = TRUE)`** ——
  读目录包装：直接读一条龙产物的标准文件名（`network/candidate_ranking.csv`、
  `metfrag/metfrag_candidates.csv`、`annotation_master.csv`、
  `multilib/<set>/{ms2_annotations_with_qc, ms1_annotations_with_confidence}.csv`、
  suspect 表），默认把四张表写到 `<run_dir>/merged/`（`candidates_ranked.csv` /
  `final_annotation.csv` / `merge_conflicts.csv` / `merge_summary.csv`）；
  `write_files = FALSE` 可只算不落盘。
- **新测试** `tests/testthat/test-annotation-rank.R`（102 断言）。

### 排序口径

- **主键 = 识别来源的置信档**，不是化合物所属的库：`L1`（in-house 命中）> `L2a`
  （online 命中）> `L2b`（预测谱库命中）> **`L2c`（仅 MetFrag 预测，无任何库命中）** >
  `L3a` > `L3b`（网络传播身份）/ `L3c`（同系物传播身份）。
  **`L2c` 是本次在本函数内部新引入的档位**（落在 `L2b` 与 `L3a` 之间），**不改动**包内
  置信度体系 `.nts_confidence_priority`（仍为 10 档，无 `L2c` / `L2d`）；in-house 化合物
  若只被 MetFrag 预测命中（库并未命中该特征），档位是 `L2c` 而不是 `L1`。
- **叠加 = 独立证据族有限晋升**，不做无界投票：候选被 **≥ 2 个独立证据族**支持时，
  `decision_label` 末尾追加 **`"*"`**（**同档、绝不越过档位**）。证据族为
  **F1 实测参考谱**（`inhouse` / `online`）、**F2 预测**（`insilico` 谱库 / MetFrag）、
  **F4 关系推断**（同系物 / 特征网络传播身份，或候选表 `ev_tier` 1–2）。
  **`ev_tier == 0` 不计** —— 它只表示「该候选本身就是一条 L1 命中」，不是独立佐证。
  MetFrag **Top-1 且 `score_total > 0`** 才参与晋升（分 0 = 库内无匹配碎片、无区分度）。
- **人工复核四条判据**（`merge_conflicts.csv`）：`mf_vs_lib_differs`（实测谱库首选与
  MetFrag Top-1 不同化合物）、`mf_vs_insilico_differs`（仅预测库时的同款对照）、
  `tied_alternatives`（同一最优档并列多个化合物）、`net_ident_differs`（网络传播身份
  ≠ 最优注释，且最优档不是实测谱档）；`conflict_status = none / resolved_by_lib`
  （最终档 L1/L2a，实测谱已定案）/ `open`（需人工看）。

### 兼容性

- **纯新增**：不动任何打分 / 过滤 / 置信度 / 网络算法，`@db_annotations` /
  `@ms1_annotations` 与导出主表、网络表**逐值不变**（0.1.31 → 0.1.32 可对旧结果做
  逐字节回归）。新函数**只读**：不写工程、不写 `ckpt_*.rds`、不写 `@network_results`。
- **导出集 202 → 204**（新增 2 个 `export()`）；est 交付脚本 `pipeline_est.R` 第 10.5 段
  改调包内 `merge_annotation_sources()`，**删除**独立脚本 `merge_rank_stage.R`。

### 补记

- 版本串（`DESCRIPTION` / `NEWS` / `AGENTS`）与内置演示脚本
  `inst/extdata/demo/demo_pfas_pipeline.R` 的 `required_version` 同步至 0.1.32。

# NTSminer 0.1.31

## 内置「三库并集去重」suspect 清单 + MetFrag 候选库默认接入

> **本版主题**：把**两个网络传播所需**、**MetFrag 候选库所需**的「三库并集去重」
> 化合物清单**固化进包内**，给出单一权威副本，免去每个脚本各自从三套库现拼。

### 新增

- **内置表 `inst/extdata/tables/NTS_PFAS_suspect_3lib_dedup.csv`（2855 行）**：
  三套谱库 `metadata` 的并集（只留 SMILES 非空），按
  `toupper(trimws(InChIKey))` **去重**，同键取「证据档更高」者
  （**inhouse > online > prediction**）；空 InChIKey 的行全部保留。
  分布 **inhouse 73 + online 103 + prediction 2679**（原始 3074 行）。
  列固定：`ID | Name | SMILES | Formula | ExactMass | InChIKey | db_type | Polarity`；
  `ID = "<集合名>_<库内 Lab_ID>"`（`inhouse_` / `online_` / `insilico_`），与
  `@db_annotations$ref_id` 同构，可回指某条库命中；`ExactMass` = 中性单同位素质量。
  安装后位于 `<library>/NTSminer/extdata/tables/NTS_PFAS_suspect_3lib_dedup.csv`，
  取值用 `bundled_nts_file("extdata/tables/NTS_PFAS_suspect_3lib_dedup.csv")`。
  口径与 `pipeline_pfas_full.R` 第 8.5 段的传播结构库**完全一致**。
  **注意**：此表只有 SMILES/质量、**无指纹**，不能直接当结构库用——要喂
  `build_structural_library(smiles_col = "SMILES", id_col = "ID")` 现场建指纹库。
- **生成脚本 `scripts/build_pfas_suspect_3lib_dedup.R`**：纯读 `.rds`、不依赖
  rcdk/JVM；写盘前硬校验「重复 InChIKey / 重复 ID / 空 SMILES / 非法质量」，
  任一非 0 即中止。

### 变化

- **MetFrag 适配层（`tests/metfrag_adapter/metfrag_adapter.R`）默认候选库**：
  `CONFIG$localdb_source` 由写死的 `NTS_Online_PFAS_OECD.rds`（159 条）改为
  **`NULL` → 自动解析包内上表（2855 条）**。解析顺序：`find.package("NTSminer")`
  → `<lib>/extdata/...`（已安装）/ `<pkg>/inst/extdata/...`（`pkgload::load_all`）
  → 再回退 `CONFIG$nts_src`；用 `find.package` 而非 `requireNamespace`，
  **不加载 NTSminer 命名空间**。仍可显式覆盖为任意 `.rds` / suspect CSV|XLSX。
  **只改默认值，不新增导出 API。**

### 兼容性

- 纯**新增内置数据 + 改一处脚本默认值**：**不动**任何打分 / 过滤 / 置信度 / 网络算法，
  `@db_annotations` / `@ms1_annotations` / 导出主表与网络表**逐值不变**，
  可对 0.1.30 结果做逐字节回归。
- **导出函数集不变**（`NAMESPACE` 无改动，仍 202 个导出）。
- 网络管线（`pipeline_pfas_full.R` / 2024SZ 交付脚本）**保持现状**——仍从加载的三套库
  逐段现拼结构库，本表只是包内权威副本；两者口径一致、可互相校验。

### 验证

- 生成脚本往返读校验：行数 **2855**、唯一 `ID`/`InChIKey` 各 2855、空 SMILES 0、
  非法质量 0；`build_structural_library()`（fread）与 `build_metfrag_localdb()`
  （`read.csv(fileEncoding="UTF-8-BOM")`）两条读取器均实测读通 2855 行。
- 版本串（`DESCRIPTION` / `NEWS` / `AGENTS`）与内置演示脚本 `required_version` 同步至 0.1.31。


# NTSminer 0.1.30

## 镜像图「两侧选错谱」修复：参考谱按命中 CE 取 + 实测谱取获胜子谱

> **本版主题**：核查「MS2 镜像图里实测与库的碎片差太多」时定位到的两个缺陷——
> 一个在**读侧**（`plot_ms2_mirror()`），一个在**写侧**
> （`annotate_by_spectral_match()`）。二者都让镜像图画出的**不是**打分时真正比较的
> 那两张谱。诊断与量化证据见
> `docs/analysis/NTSminer_镜像图两侧选谱缺陷_0.1.30.md`。

### 缺陷与修复

- **D1（读侧，影响 72/72 命中）参考谱跨 CE 合并**：`.nts_get_ref_ms2()` 把
  `@spectra_long` 里该化合物**所有**碰撞能量（CE10 + CE20 + CE40 …）round-4 后按
  `max(intensity)` 合并成一条「共识谱」当作参考半谱。于是副标题写着 `lib_CE CE10`，
  下半谱却是三档 CE 的叠加——PFOA（FT2224）实测 3 峰对上的 `CE10` 档（库中亦 3–6 峰）
  被画成 30 峰的全 CE 叠加，视觉上「差太多」。修复：`.nts_ref_levels_from_flat()`
  新增 `nce` 过滤（`.nts_get_ref_ms2()` 透传 `lib_CE`），参考半谱**只取命中的那一档
  CE**；该档缺失或旧工程无 `lib_CE` 时**回退**到全 CE 共识，行为不变。legacy
  `databaseClass` 路径同样按 `lib_CE` 收窄（仅当嵌套列表以 CE 命名时）。
- **D2（写侧 + 读侧，影响 57/72 命中）实测谱取错子谱**：一个特征常聚到多个 MS2 扫描
  （多文件 / 多 CE），`annotate_by_spectral_match()` 在所有 `(库 CE × 实测子谱)` 组合里
  选**复合分最高**的获胜子谱打分，却**没有把获胜子谱的索引写进 `@db_annotations`**；
  `plot_ms2_mirror()` 于是回退到 `idx = 1`，画的是**第一个**子谱，而非参与打分的那张。
  修复分两步：
  - **写侧**：`annotate_by_spectral_match()` 记录 `best_ms2_index`（获胜实验谱在
    `@ms2data` 子谱列表中的位置）与 `n_ms2_subset`（该特征提供的子谱数），
    写入 `@db_annotations`；
  - **读侧**：`.nts_get_best_exp_ms2()` 优先用显式 `best_ms2_index`，没有时按
    `exp_CE` 标签（如 `CE10#1835@90ppb2`）匹配子谱，都没有才回退 `#1`。标签匹配让
    **旧工程（0.1.27 及更早，无 `best_ms2_index` 列）也能选对**获胜子谱。

### 镜像图：匹配碎片标注 m/z（绘图默认开启）

- `plot_ms2_mirror()` 新增 `label_matched_mz = TRUE` / `label_mz_tol = 0.02` /
  `label_digits = 2`：**两半都画出来**的匹配碎片离子，各自在峰尖旁标出自己的 m/z
  （实测半红字、库半蓝字，均竖排）。匹配**就在画出的峰上做**（扣噪后、容差 0.02 Da，
  与打分口径一致），所以标注的峰**一定在图上**，不会出现「标了却看不到的峰」。
- 相邻太近（不足一个字宽）的匹配峰会**自动错开到第二层竖直位置**，数字互不压叠；
  因此图的上/下留白会按需加大。`label_matched_mz = FALSE` 回到原样（y 轴 ±105、
  无标注），旧行为零改变。
- 纯**显示层**改动：不触碰任何打分 / 过滤 / 置信度 / 导出数值。`export_mirror_report()`
  同名三参数透传。

### 兼容与回归

- `@db_annotations` **新增两列**（`best_ms2_index`、`n_ms2_subset`），**不删列、不改
  数值**：对同一工程重跑 `annotate_by_spectral_match()`，新旧输出在全部 27 个旧列上
  **逐位一致**（实测 max |Δ| = 0），仅多出这两列。
- 镜像图**视觉变化、数值不变**：`score_total` 等打分列不受影响；改的是图上画哪两张谱
  （并新增匹配碎片 m/z 标注）。无 `lib_CE`/`exp_CE`/`best_ms2_index` 的旧工程自动回退旧行为。
- **副标题**由 `exp_CE CE10 | lib_CE CE10` 改为 `exp CE10#1835@90ppb2 | lib CE10`，
  直接标出双方用的**具体**实验扫描与库 CE 档。
- 本版**未改**任何打分、去重、置信度或网络逻辑。

# NTSminer 0.1.29

## 全流程逐函数核查修复：2 个 P0 + 12 个 P1

> **本版主题**：一次覆盖**全包 48 个文件 / 670 个函数**的逐函数核查（按数据处理流程
> 逐段核对签名与数据契约）产出的缺陷清单，本版收口 **2 个 P0 + 12 个 P1**（14 个
> 条目，落在 **12 个文件**）。核查过程与逐条证据见
> `docs/analysis/NTSminer_全流程函数核查_0.1.28.md`（含 §8 修复记录）。

### 修复：2 个 P0（静默数据错误）

- **P0-1 `.name_from_ref` 作用域缺失崩溃**：函数体引用了未在其作用域内绑定的列，
  在结构传播/惰性求值路径下会抛「object not found」。改为显式取列，消除符号逃逸。
- **P0-2 GC 跨文件 `compound_id` 重号致强度错位**：GC 数据合并多文件时按文件名计数，
  `compound_id` 在不同文件间重复，导致对齐强度矩阵时把**不同化合物**叠加到同一行。
  改为全局唯一编号（含文件维度），强度不再错位。

### 修复：12 个 P1（口径/边界/静默丢行）

- **P1-1** 写出/导出路径的 `data.table` 右连接残留 **`i.*` 前缀列**污染输出表 ——
  在安全点清理，导出表列名回归预期。
- **P1-2** 同位素比值校验：**无法计算**的比值（`NA`）此前被当作 `FALSE` 丢弃，
  与日志「跳过计数」语义矛盾。改为 **`NA` 视为接受并计入 `n_ratio_skipped`**，
  只有**显式计算失败**（`FALSE`）才丢弃。
- **P1-3 KMD 默认窗口不可达**：`calculate_kmd` 的 KMD ∈ [-0.5, 0.5]，而默认
  筛选窗写成 `c(0.85, 1)`（[0,1) 约定），**恒为空**。改为镜像为 `c(-0.15, 0)`，
  与 `c(0, 0.15)` 合并覆盖 |KMD| ≤ 0.15，对齐 `calculate_kmd` 自身默认
  `md_range = c(-0.15, 0.15)`（**行为变化**）。
- **P1-4 同位素理论峰间距未按电荷缩放**：Rdisop 返回**中性**单同位素质量，多电荷
  离子峰间距应按 ~1/|z| 缩放。新增 `charge` 参数（默认 1 为 no-op，保持中性质量
  行为不变）；加合物偏移仍由调用方先对齐 M+0。单峰 `cosine = 1` 分支不动（测试钉死）。
- **P1-5 DIA 前体 EIC 误用碎片阈值**：`min_frag_intensity` 是**碎片**顶点阈值，
  却被用来闸前体 EIC（量纲不同），导致本可命中的前体被砍。改为仅要求前体 EIC
  **≥3 个扫描点有信号**，真正的甄别交给后续共流出检验（**行为变化**）。
- **P1-7 置信度库列名判空残缺**：`db` 表列名探测漏了 `feature_id` 分支，某些库
  走到取 `name` 列时取空。补齐 `name → feature_id` 回退。
- **P1-8 内标 RT 分段映射越界**：`seg_is[seg]` 用**段号**去索引**每个内标**的向量，
  段数与内标数不等时越界/错配。改为 `match(seg, seg_is)`，无内标的段返回 `NA`。
- **P1-9 / P1-9b RT 与 logIE 模型评估行错位**：`predict_rt()` / `predict_logIE()`
  会**静默丢行**（NA/空 SMILES，经 `compute_fingerprints()` 的 `valid_mask`），
  而评估函数按**原长**比对，令观测与预测**错位**。改为按同一 `smi_ok` 掩码重展开，
  长度不符时显式 `warning` 不再静默。
- **P1-10 描述符零方差列剔除误伤单行批次**：非 NA 值 **< 2 个**的列**无可估方差**，
  此前被当零方差删掉，单行批次会丢掉**所有**描述符列。改为仅当 `length(xc) ≥ 2`
  且 `var < 1e-12` 才删。
- **P1-11 `.predict_rie_df` 丢弃全部 logIE 预测**：同 P1-9 的行数错配使余下用
  mock 兜底值填 RIE。改为按 `smi_ok` 重展开，长度不符时 `warning` 并回退。

### 未变更（经复核判定为设计如此）

- **P1-6（原列 P1）：`batch_annotate()` 按 `ref_id` 一对一去重**。复核调用方
  `06_scoring_engine.R`（含 `@ms2_score_grid` 的 `is_winner` 同步注释）确认这是
  **有意且已文档化**的行为，改动会改变已验证基线（in-house 70 → 93）。**保持现状**，
  降级为「设计如此 / 待产品决策」。

### 说明

- 本版全部改动为**盲改 + 静态核对**（核查环境无 R），**未运行** `R CMD check` /
  单元测试 / 回归实跑；使用外部工具的批次（RT/logIE/内标/KMD）建议在本地装包后
  复跑一次验证。
- P2（约 40 条批次）**未处理**。

# NTSminer 0.1.28

## 内置数据可用化：常用 PFAS 谱库 + 无参演示数据 + 统一查找链

> **本版主题**：让**新装的 NTSminer 不依赖任何外部文件即可跑通全流程**。此前包内的
> 表 / 模型 / 谱库散落在 `inst/extdata` 与被 `.Rbuildignore` 排除的顶层 `data/`，
> 且代码用 `system.file("data", ...)` 去取 —— **装上包后恒取不到**，默认参数指向的
> 文件实际全部缺失。本版把「包内自带文件」收敛成**单一查找链 + 固定目录约定**，
> 并把常用 PFAS 谱库与一套合成演示数据真正放进包里。

### 新增：统一的内置文件查找链与内置谱库

- **`.nts_extdata()`（内部）**：所有「默认指向包内文件」的解析点走**同一入口**，
  依次尝试 `<pkg>/extdata/<rel>`（已安装）、`<pkg-root>/inst/extdata/<rel>`
  （源码树 / `pkgload::load_all()`）、`<cwd>/extdata/<rel>`，返回**第一个存在**的
  规范化绝对路径；支持**新旧布局并存**（新 `tables/…` 优先，兼容旧根路径）。
  找不到时 `mustWork = FALSE` 返回 `NA_character_`（调用方优雅降级），
  `TRUE` 则给出**含期望路径的清晰报错**。
- **导出 `bundled_nts_file(rel, mustWork = TRUE)`**：查找链的用户可见包装，
  便于拿包内表的磁盘路径，如
  `bundled_nts_file("extdata/tables/diagnostic_fragments.csv")`。
- **目录约定（相对 `inst/extdata/`）**：`libraries/` 常用谱库、`tables/` 专家表、
  `models/` logIE / RT 模型、`demo/` 演示数据。
- **导出 `nts_builtin_library_paths(which)` / `load_builtin_pfas_library(which, combine)` /
  `merge_nts_libraries(dbs)`**：随包内置**三套常用 PFAS 谱库** —— in-house 标品库
  （`inhouse`，置信 L1）、OECD online 库（`online`，L2a）、OECD in-silico 库
  （`prediction`，L2b）。`load_builtin_pfas_library()` 默认把三库**合并成一个
  `NTSDatabase`**（逐行保留各自 `db_type`），一次 `batch_annotate()` 即可同时触达
  L1/L2a/L2b；单库请求原样返回（不加 `Lab_ID` 前缀）。**超大社区库
  （MoNA-unified ~15 MB、spectraverse ~600 MB）不内置**。
- **`predict_rie()` 的 in-house logIE 模型**与 **`default_diagnostic_fragments_csv()`**
  改走查找链，模型现落在 `extdata/models/logIE/`、专家表落在
  `extdata/tables/diagnostic_fragments.csv`（旧根路径仍兼容）。

### 新增：无参演示数据（全流程零外部文件可跑）

- **`inst/extdata/demo/`** 自带一套**合成** LC-MS/MS 数据（负离子、PFAS 类似物）：
  `ms1/`（`Sample_1..3` + `Blank_1..2`，5 个 centroid mzML）、`ms2/`
  （`Sample_1..3`，3 个 DDA mzML），合计约 1.7 MB。化合物品名、保留时间、母离子
  m/z 与 MS2 碎片**全部取自包内 in-house 库**——只有**原始信号**是模拟的，因此演示
  确实是拿**随包发布的 L1 参考数据**做注释。数据由包根 `tools/make_demo_data.R`
  **确定性生成**（无 RNG、无时间戳，重跑逐字节一致）。
- **`inst/extdata/demo/demo_pfas_pipeline.R`**：可直接
  `source(system.file("extdata", "demo", "demo_pfas_pipeline.R", package = "NTSminer"))`
  的**无参**端到端脚本。它是生产脚本 `tests/rawdata_pfas/pipeline_pfas_full.R` 的
  **内置数据版**：段序、函数调用、导出表完全一致，只有顶部 CONFIG 指向包内数据
  （`MS1_DIR` / `MS2_DIR` / `db_sets` 都用 `system.file()` 取路径），因此**演示脚本的
  写法就是用户自己脚本的写法** —— 换成自己的数据只需改 CONFIG。**`DEMO.md`** 记录
  数据说明、段序表、期望结果与参数注意事项（`peakwidth` 单位是**秒**；诊断碎片表是
  PFAS 专属规则、需显式传入）。
- **演示脚本带前置自检**：脚本第 0 段核对**本脚本真正用到的每一个函数是否存在**
  （`load_nts()` / `have_api()` / `required_fns`），而不只看 `packageVersion()`。
  `packageVersion()` 读的是**磁盘上**（或源码树）的 `DESCRIPTION`，可能比当前
  **已挂载**的命名空间新 —— 例如源码树改完没重跑 `load_all()`、或 `R CMD INSTALL`
  后**没重启 R**、或 `.libPaths()` 里有个更旧的副本抢先命中。此时旧版会以
  `没有"xxx"这个函数` 这种**看不出原因**的形式失败；现在改为**明确列出缺少哪些函数**、
  并提示把 `DEV_SRC` 指向源码目录（或设 `NTSMINER_DEV_SRC` 环境变量）回退源码加载。
  CONFIG 第 4 块另有一处内置数据自检：`system.file()` 取不到演示数据/谱库时直接报出
  `find.package()` 命中路径与重装命令。
  （说明：本版早期曾导出 `nts_demo_ms1_dir()` / `nts_demo_ms2_dir()` /
  `nts_demo_project()` 三个「演示专用」取路径函数，现已**移除**——演示改为走正常
  的数据处理流程，不引入只有演示才有的 API。`bundled_nts_file()` /
  `nts_builtin_library_paths()` / `load_builtin_pfas_library()` / `merge_nts_libraries()`
  保留导出。）

- **实测结果**（单线程 ~50 s，本条为**全流程**跑完的数字）：38 个 MS1 特征 →
  空白扣除后 35 → 22 个同位素标记 → 33 个特征落在 CF2 KMD 序列 → MS2 匹配 11 个
  → 净化后 9 张干净谱 → 无库 CF2 阶梯成立 5 个特征（9 个有子类提示）；
  逐库注释：in-house **12 行 / 11 特征（全 L1）**、online 5 行 / 4 特征（L2a）、
  in-silico 66 行 / 11 特征（L2b）、合并 **83 行 / 11 特征**；合并网络 38 节点 /
  20 边，结构指纹库现场构建（`structural_db.rds`）并出 4 张交互 HTML；
  `annotation_master.csv` 35 行。11 个加标化合物全部正确且为 **L1**；in-house 第 12 行
  是 PFHxPA 对 PFHxS 特征的**等重近邻误命中**（m/z 差 0.3 ppm，`score_total` ≈ 0.67，
  因打分为加性而非 RT 硬闸），特意保留以展示真实边界情形。
  > 一致性提示：`filter_blanks()` 裁剪 `feature_table`（38 → 35）后，脚本第 4 段会把
  > `ms1_annotations` **同步裁到同一批保留特征**，故 `ms1_annotations_with_confidence.csv`
  > 与主表都是 **35 行**（本版修复，见下一条）。

### 修复：空白过滤后 ms1 视图未同步裁剪（脚本）+ 去掉 `@network_results` 遗留 stub（包内）

- **脚本第 4 段「空白过滤后同步 `ms1_annotations`」此前是空操作**：
  `get_ms1_annotations(x)`（`step = NULL`）返回的是**扁平 `data.frame`**（一行一特征、
  键 `name`），而旧代码写成 `lapply(ms1, function(t) if (is.data.frame(t) ...))` —— 实际在
  **遍历列**，于是什么都没裁掉，空白过滤后 `ms1_annotations` 仍带着被扣掉的特征行
  （演示数据：特征表 35 行，`ms1_annotations` 却 38 行）。现改为**直接按 `name` 过滤扁平表**
  再 `set_ms1_annotations()` 写回。同步修改 `tests/rawdata_pfas/pipeline_pfas_full.R`、
  `tests/rawdata_pfas/standard_pipeline_MS1_MS2_to_annotation.R` 与内置演示脚本
  `inst/extdata/demo/demo_pfas_pipeline.R`（演示脚本主体与生产脚本**逐字节一致**）。
- **`@network_results` 不再被注释步骤写入**：`annotate_by_spectral_match()` 历史上会往
  `@network_results` 塞一个 **legacy stub**（`nodes` = 每特征最佳行、`edges` = 空表），
  是分子网络模块之前的遗留物，与真正的 MS2 相似度网络（`build_molecular_network()` /
  `get_network_data()`）无关 —— 只会让挖掘后的工程莫名带上一份非网络内容。现**移除该写入**，
  注释过程只写 `@db_annotations` / `@ms2_annotations`，`@network_results` 恒为空
  （除非调用方显式 `set_network_results()`）。**槽位与 `get_`/`set_network_results()`
  访问器保留**，以免破坏已保存工程的 `readRDS` / `validObject`。同步更新
  `tests/testthat/test-annotate-per-ce-grid.R`（改为断言该槽为空）。

### 修复：`calculate_exact_mass()` 不再硬依赖可选的 `Rdisop`

- 该函数被 **KMD 基准质量、ISF 中性丢失、同系物 CF2 步长**共用，却只放在
  `Suggests` 的 `Rdisop` 里实现：**未装 `Rdisop` 时全部返回 `NA`**，直接导致
  `flag_isf()` 抛 `neutral_losses must contain at least one finite value`、
  `calculate_kmd()` 无法建立序列。现改为**双引擎**：
  1. **先试 `Rdisop`**——已装则数值与历史版本**逐位一致**；
  2. 未装则用**包内置单同位素质量表**（`.NTS_MONOISOTOPIC_MASS`，37 个元素）配
     手写 Hill 式解析器（`.nts_formula_counts()`，支持 `"C3HF5O2"`、
     `"(CF2)3"` 等含括号+系数的写法）计算，**结果与 `Rdisop` 一致到 <1e-6 Da**；
     首次替代时**一次性 `message()`** 说明并提示安装 `Rdisop` 可取逐位一致值。
- 影响：**已装 `Rdisop` 的环境数值不变**（因此可用旧基线回归）；最小安装下
  KMD / ISF / 同系物步长**由「整段不可用」变为「可用」**。

### 其他

- **Shiny 工作台计入构建**：`.Rbuildignore` 移除 `^inst/shiny_app`。此前 `app.R`
  被排除在构建外，`launch_ntsminer()` 装包后必然报
  `Could not find shiny_app directory`（并误导用户「重新安装」）。`app.R` 只读
  用户上传文件、无其它外部依赖，纳入构建即可修复。
- **文档**：vignette「数据准备 / Phase 1 / Phase 3」由失效路径
  （`../inst/extdata/ms1|ms2`、`readRDS("../inst/extdata/database/…")`、
  `diagnostic_csv_path = "../inst/extdata/diagnostic_fragments.csv"`）改为
  `system.file("extdata", "demo", "ms1"|"ms2", package = "NTSminer")` /
  `load_builtin_pfas_library()` / `default_diagnostic_fragments_csv()`，并补
  `peakwidth` 单位与演示数据说明。
- **同步**：`NAMESPACE` + 手写 `man/{bundled_nts_file,nts_builtin_libraries}.Rd`、
  `man/calculate_exact_mass.Rd`；`R/05_ms2_series.R`、`R/08_risk_assess.R` 中过期的
  `data/…` 路径注释更正。**未改动任何过滤 / 打分 / 置信度算法与 `@db_annotations` 口径**。

# NTSminer 0.1.27

## 镜像谱图（MS2 mirror plot）修复：适配扁平库 + 噪声线对齐打分口径 + 得分拆解

> **行为变更**：`plot_ms2_mirror()` / `export_mirror_report()` 现可接受**扁平 `NTSDatabase`**
> 参考库（`@spectra_long`，即 `lc_libraries/` 下全部库的格式），不再只支持老式 `databaseClass`；
> 噪声下限 `min_rel_intensity` **默认 `2.5`（基峰的 2.5%，与打分引擎 `int_thresh = 0.025`
> 去噪口径一致）**，出图前丢弃低于该阈值的峰；副标题新增 **`score_total` 与 `MS2` 的得分拆解**。

- **修复**：`.nts_get_ref_ms2()` 原先硬性要求 `database` 是 `databaseClass`（嵌套
  `@spectra.data$Spectra.negative[[Lab_ID]]`），对**实际在用的扁平库**（in-house / online /
  prediction / merged / MoNA / gnps…）直接报
  `database must be a databaseClass for mirror plotting`。现改为**双栈分派**：
  `NTSDatabase`（扁平——从 `@spectra_long` 按 `Lab_ID` + `Polarity` 取行、按 `NCE` 分组还原成
  `list(NCE -> data.frame(mz, intensity))`）与 `databaseClass`（老嵌套）**都支持**；列名做
  `Lab_ID`/`lab_id`、`frag_mz`/`mz`、`Polarity`/`polarity`、`NCE`/`nce`/`CE` **容错**，缺 `NCE`
  时折叠为单一 `Merged` 层；`ion_mode = "auto"` 仍先负后正。参考谱口径
  （`.nts_lib_consensus_matrix`：m/z 4 位分箱取各 NCE 最大值）**未变**。
- **噪声下限 `min_rel_intensity`（默认 `2.5`）**：两半谱各自归一到基峰 = 100 后，丢弃
  `|intensity| < min_rel_intensity` 的峰，去掉低强度「噪声地毯」；传 `0` 保留全部（旧行为）。
  **默认值取 `2.5` 是为了与打分器一致**——打分引擎 `.calculate_clean_ndp()` 以
  `int_thresh = 0.025`（基峰 2.5%）去噪，故镜像图现在只画*打分器真正参与比较的峰*，
  不会再出现「镜像看着不像、分数却很高」的口径错位。
- **副标题新增得分拆解**（内部助手 `.nts_score_breakdown_lines()`）：第一行
  `score_total = m/z + RT + MS2`（按 `0.25 / 0.25 / 0.50` 权重加权后的三项），第二行
  `MS2 = ndp + nl + md + diag`（复合四项，附 `raw ndp`、匹配峰数与 `score_diagnostic` 原始诊断覆盖率）。
  其中复合项的 `diag` 取残差 `max(0, MS2 − ndp − nl − md)`（扁平库存的 `score_diagnostic` 是**原始覆盖率**、
  未进复合式，故额外单列显示以便对账）。如此一眼可见
  `score_total` 并非纯谱图相似度——MS1（母离子 m/z + RT）两项常占约一半。列缺失时显示 `NA`，
  自动判别新旧两种加权口径（扁平库 `0.25/0.25/0.5`，老 `databaseClass` 存的是已加权值）。
- **影响**：只影响绘图函数与其新增参数，**不动 `@db_annotations` / 打分 / 置信度**，不写 `prj`。
- **同步**：`man/{plot_ms2_mirror,export_mirror_report}.Rd`；PFAS 一体化脚本
  `pipeline_pfas_full.R` 与 2024SZ 副本 `pipeline_pfas_2024SZ.R` 的 `required_version`
  对齐到 `0.1.27`。

# NTSminer 0.1.26

## 默认参数变更：`min_frag_freq` 由 `0.8` 改为 `0.5`

> **行为变更（默认值）**：`generate_network_data()` / `get_network_data()` /
> `plot_nts_network()`，及其内部编排 `.nts_build_network_edges()` /
> `.nts_ms2_layer_edges_dt()`，形参 **`min_frag_freq` 的默认由 `0.8` 改为 `0.5`**。
> 仍为**只读派生**，**不写 `prj` / `@network_results`**；版本仅此一处语义变化，无 API 增删。

- **背景**：`min_frag_freq` 是 `@clean_ms2` 共识峰的 **`frequency` 下限（floor，取值 0–1，
  非检出次数）**，在建 MS2 特征网络（DFI / NL / MD）前逐峰过滤。真实 PFAS 数据实测 `0.8` 过严——
  966 个特征里仅 633 个保留 ≥2 个峰、69 个库种子仅 38 个够成边；`0.5` 时成边种子明显增多，
  同时仍滤掉只出现在 <1/5 扫描里的杂散离子。（检出次数闸是另一处 `purify_min_occurrence`，未改。）
- **影响**：以默认参数直调 `get_network_data()` / `generate_network_data()` /
  `plot_nts_network()`（如 Shiny App / `plot_molecular_network` 路径）时 **MS2 边数量与构成会变化
  （边增多）**。**需保留旧默认者请显式传 `min_frag_freq = 0.8`**；`NULL` / 非有限值仍表示不过滤。
- **三档对照**（真实 PFAS 工程，只读，见 `tests/rawdata_pfas/run_manual/network_from_prj/compare_nf/`）：
  成边种子 0.8 / 0.5 / 0.1 = DFI 29/54/62、NL 7/29/49、MD 11/39/61、合并 61/67/69；
  边数 0.8→0.5 放大 3–7×，0.5→0.1 再翻 2–3× 但收益饱和，故取 **0.5**。
- **同步**：`man/{generate_network_data,get_network_data,plot_molecular_network}.Rd` 默认值文本；
  PFAS 流程脚本 `standard_pipeline_MS1_MS2_to_annotation.R`（`network_min_frag_freq`）与
  `pipeline_from_project.R`（`net_min_frag_freq`）取 `0.5` 与包内默认对齐。

# NTSminer 0.1.25

## 分子网络（特征网络，路径 A）S1–S3 加固：词汇统一 / 边权归一 / 规模护栏 / 同系物步长 / CF2 阶梯加固

> 全部为**只读派生**改动：`get_network_data()` / `generate_network_data()` 依旧**只返回**
> `list(nodes, edges)`，**绝不写 `@network_results` 或任何 `@` 槽**；新增内容均为可选列或
> 带旧默认的可选参数（默认 no-op），保证既有调用与导出口径不破坏。
> 设计见 `docs/analysis/NTSminer_分子网络设计方案_v0.1.24.md`。

- **S1-a 边词汇统一（P1-5）**：`@` 无关，`ms2_edge_types` / `ms1_edge_types` 的别名解析
  现**大小写不敏感**。MS2 层别名：`dfi`→`fragment`、`nl`→`neutral_loss`、
  `md`/`massdiff`/`mass_diff`→`mass_difference`；MS1 层别名：`kmd`/`hom`/`homolog`/`series`→
  `homologue`、`trans`/`transform`→`transformation`。
- **S1-a `Transformation` 独立成类 + 新增 `type="transformation"` 切片**：
  `edge_category` 不再把 `Transformation` 混并入 `Homologue`，而是单列 `"Transformation"`；
  `get_network_data(type=)` / `.nts_network_edges_filter_by_type()` 新增 `"transformation"` 档位。
  ⚠️ **行为变更**：`type = "homologue"` 现**只**返回同系物边（此前会把转变产物边一并返回）。
  需要旧并集者请显式取 `type="all"` 或分别取 `"homologue"` + `"transformation"`。
- **S1-b `edge_score` 归一化（P1-4）**：MS2 计数边（DFI/NL/MD）的 `edge_score` 由过去的**恒为
  1.0** 改为**归一化权重** `n_shared / min(peaks_a, peaks_b)` ∈ (0, 1]，越大代表共享证据越强
  （「2 个峰里共享 2 个」> 「40 个峰里共享 2 个」）。`min_ms2_score` 默认由 `0.6` 改为
  **`NULL`（不过滤）**，以保持旧默认行为；需要按分数筛边时显式给数值（如 `0.5`）。
  ⚠️ 依赖 `edge_score` 做加权的下游（`annotate_by_dual_network` 的传播比较器）排序会更细，
  但**无硬阈值**，不改变可行解集合。
- **S2-a 规模护栏（P0-3）**：`get_network_data()` / `generate_network_data()` 新增 `top_k`
  （每个节点最多保留的边数，按 `edge_score` 从强到弱）、`max_edges`（全图边数上限，按分数截断）、
  `min_degree`（丢弃度低于此的节点）。**默认 `Inf` / `0L` 即 no-op**，旧调用结果不变。
- **S2-b 同系物边补齐步长与标签 + 精确 CF2 常量（P1-1）**：MS1 `Homologue` 边新增
  `mass_difference`（相邻同系物 m/z 差）、`n_steps`（`round(Δm/z / CF2)`）、`transformation`
  （`+nCF2` 带符号标签），`label` 由泛化的 `KMD_Homologue` 改为 `+nCF2`（超出一个 CF2 阶梯
  容差的记为 `""`）。同时把 `build_homologue_network()` 默认 `target_mz_diff` 由舍入值
  `49.9968` 统一为**精确常量 `49.9968064`**（`.nts_pfas_cf2_da`），`.mn_precursor_delta_rule_table()`
  的 CF2 项同步改为引用该常量。
- **S2-c `md_min_da` 参数化（P1-2）**：MD（质量差）边的谱内碎片间距下限由**硬编码 1.0 Da**
  改为参数 `md_min_da`（默认仍 `1.0`，等价旧行为）。
- **S3-a CF2 阶梯入边加固 `cf2_ladder_boost`（默认开）**：两端特征都带 MS2 CF2 阶梯
  （节点表 `pfas_ladder_len >= 3`，来自 0.1.23 的 `detect_pfas_fragment_series`）的 `Homologue`
  边，标记 `ms2_ladder_confirmed = TRUE` 并赋 `edge_score = 1`。**只加固已有 MS1 同系物边，
  不新增特征级边**（与 0.1.23「CF2 阶梯只作证据、不把碎片当节点」的裁定一致）；节点表无
  `pfas_ladder_len` 时为 no-op。
- **S3-b 阶梯可视化（P2-3）**：`plot_molecular_network()` 节点 tooltip 增加 **PFAS ladder**
  行，展示 `pfas_ladder_len >= 3` 的 CF2 阶梯长度（含 `pfas_subclass_hint` 子类提示）。

---

# NTSminer 0.1.24

## 改名：链长两列 `pfas_chain_lo` / `pfas_chain_hi` → `pfas_chain_from_series` / `pfas_chain_from_ladder`

- **`@ms2_annotations` 两列改名**（0.1.23 新增的 9 个 `pfas_*` 列里的最后两列）：
  `pfas_chain_lo` → **`pfas_chain_from_series`**，`pfas_chain_hi` → **`pfas_chain_from_ladder`**。
  **数值本身不变**，只是名字与含义对齐；其余 7 列、`pfas_series_score`、`is_hit` 判据、
  `[MS2 Verified]` 标记均不受影响。
- **为什么要改**：这两列按设计 §4.5 是**两个独立口径、互为校验**，并非「下界 / 上界」区间——
  `from_series = max(1, pfas_series_max_n)`（直接观测到的 CnF(2n+1)- 系列长度），
  `from_ladder = pfas_ladder_span + 1`（由 CF2 阶梯跨度外推）。任一都可能更大，实测出现过
  逆序（如 FT3015 `from_series = 10` / `from_ladder = 6`，即旧名下的 `lo > hi`），`lo` / `hi`
  容易被误读成区间。
- **`find_pfas_fragment_series()` 返回值的 `chain` 元素同步改名**：
  `chain$lo` / `chain$hi` → `chain$from_series` / `chain$from_ladder`（详见 `?find_pfas_fragment_series`）。
- ⚠️ **下游按旧列名取列的代码须改**（如 `m2$pfas_chain_lo`）。已有的 checkpoint 与已导出的
  `ms2_annotations_with_qc.csv` / `pfas_series_summary.csv` 里仍是旧列名——**删掉断点重跑**即可换成新名。

---

# NTSminer 0.1.23

## 无库 MS2 CF2 碎片阶梯 / PFAS 系列检测（分子网络后续 P0-1）

- **新增纯函数 `find_pfas_fragment_series()`**：对单张 MS2 谱做无库结构判定——
  ① 在排序后的 m/z 轴上求**最长 CF2 阶梯**（相邻峰相差整数倍 CF2 = 49.9968064 Da；
  `mode = "loose"` 允许缺环、单跳最多 `max_gap_units` 个 CF2，`mode = "strict"` 只认单步）；
  ② 匹配全氟烷基负离子系列 CnF(2n+1)-（CF3- 锚点 68.9957596，步长 CF2）；
  ③ 匹配内置 PFAS 诊断碎片专家表（`default_diagnostic_fragments_csv()`）；
  三者合成 **子类提示**（PFCA / PFSA / FTS）与一个连续 **证据分**。
  返回 `list(ladder, series, diagnostic, score, chain, is_hit, ...)`，无工程对象、无 I/O 副作用，专为单测设计。
- **新增 `pfas_diagnostic_rules()`**：内置化学规则表（P1-2），把诊断离子质量映射到 PFAS 子类。
  PFCA = 任一 CnF(2n+1)-；PFSA = 任一磺酸 / 氟磺酸 / 三氟甲磺酸头基离子；FTS = 磺酸头基**且**存在 CF2 阶梯
  （`type = "series"` 且 `required = TRUE`）。子类成立需「该子类至少一行命中」且「所有 `required = TRUE` 行命中」。
- **新增 `detect_pfas_fragment_series()` 并接入主流程（默认开，`detect_pfas_series = TRUE`）**：
  `annotate_by_spectral_match()` 与 `batch_annotate()` 现默认调用之，对每个特征的一致 MS2（`@clean_ms2`）
  跑检测并把设计 §4.6 的 **9 个新列**写入 `@ms2_annotations`：
  `pfas_ladder_len / pfas_ladder_mz / pfas_ladder_span / pfas_series_max_n / pfas_series_hits /
  pfas_subclass_hint / pfas_series_score / pfas_chain_lo / pfas_chain_hi`。
  **只新增列，绝不改写既有列**。
  - 联合证据门槛（设计 §8.7）：`is_hit = 阶梯成立 AND (CnF(2n+1)- 命中 >= 1 OR 子类非空)`，
    避免仅凭「巧合的 ~50 Da 等间距」误判。
  - 命中特征的 `@ms1_annotations$annotation_level` 打上 `[MS2 Verified]`，与 0.1.21 诊断碎片走同一条代码路径。
  - `@clean_ms2` 为空时为安全 no-op（仅记录日志步骤）。
- **分子网络：只做节点属性（决策 #3）**。新列经既有 `@ms2_annotations -> 网络节点` 合并自动流入节点，
  **不新增边类型**；`pfas_series_score` 不匹配节点的 `^score_` 剪列规则，稳定保留。
- **置信度结构护栏加入 `PFAS_Series` 标记（决策 #2）**：`assign_confidence_levels()` 的
  `annotation_method` 结构判定现额外识别 `PFAS_Series`。命中特征若其 `@db_annotations$annotation_method`
  **本就带结构 token**（`Dual_Network` / `Propagation` / `Network_*` / `Homologue` / `Transformation`），
  追加 `; PFAS_Series`；该操作**可证明地不改变置信档**（该行本就按结构处理）。非结构行与空行绝不改动，
  **不新增置信级别、不做提级**。
- 新增单元测试 `tests/testthat/test-pfas-fragment-series.R`（纯函数、子类语义、正负极性、工程级集成、结构护栏）。
- 设计文档：`docs/analysis/NTSminer_MS2_CF2碎片阶梯检测设计方案_v0.1.22.md`（S1–S4 已实现）。


# NTSminer 0.1.22

## 置信度分级改为「按库来源」划分（破坏性变更）

- **三个库档不再按 RT / MS2 余弦细分，直接由谱库来源（`db_type`）决定**：
  `in-house → L1`、`online → L2a`、`prediction → L2b`。之前 in-house 有无 RT 分 L1/L2、
  online 有没有过余弦门槛分 L2a/L2b、prediction 有无 RT 分 L2c/L2d，现全部合并为「一个来源一档」。
  - `assign_confidence_levels()` 的库类型分支简化为 `if inhouse L1 / online L2a / prediction L2b`，
    不再读取 `rt_error_s`、也不再检查 `score_cosine_ms2`。
  - **`L2`（in-house 无 RT）、`L2c` / `L2d`（prediction）随之退场**，方案由 **13 级缩为 10 级**：
    `L1, L2a, L2b, L3a, L3b, L3c, L3d, L3e, L4, L5`。`.nts_confidence_priority` 与
    `.nts_confidence_label` 同步收缩，下游 L3a 起的优先级整数相应下移。
  - RT 一致证据仍保留，但**只用于**网络传播档 `L3a`（propagation + RT）与 `L3b`
    （propagation − RT）的区分；`rt_tol_sec` 现仅影响该档。
  - `ms2_cosine_threshold` 参数**弃用**（保留形参，避免破坏既有调用；不再起作用）。
- **修复：`batch_annotate()` 现把库来源写进 `@db_annotations$db_type`。** 此前只有
  `annotate_by_spectral_match()` 带上了该列，而主流程入口 `batch_annotate()` 自行重建
  `db_annotations` 时漏掉，导致自建库命中无法被判为 in-house、经兜底逻辑一律按 online 计
  （in-house 标品命中被降成 L2a，且触发「no `db_type` column」警告）。现在逐行从（已盖来源戳的）
  MS1 候选表 `comp` 回填 `db_type`，in-house 命中稳定落到 **L1**。
- 只有**真正的库命中**才按库来源分档：`@db_annotations` 中 `annotation_method` 标记为**结构性**结果
  （网络传播 `Dual_Network`/`Propagation`/`Network_*`、同系物 `Homologue`、转化 `Transformation`）
  的行**不参与**库分档，保留其 `L3a`/`L3b`/`L3c`/`L3d`/`L3e`/`L4` 档，避免被附带写入的 `db_type`
  误标成 L1/L2a/L2b。（同时含 `Spectral_DB` 的混合行仍视为库命中。）
- 影响：`confidence_level` 取值与计数会变化（in-house 命中一律 L1、online 一律 L2a、
  prediction 一律 L2b），单元格标签与因子水平集都变了；依赖旧标签的下游/报告须同步。设计取舍是
  **置信度只反映「命中来自哪种库」**，不再反映单次匹配的 RT/谱图强弱（弱命中已被 `batch_annotate()`
  的 `min_score` 门槛挡在 `db_annotations` 之外）。如需回退旧行为，请锁定 0.1.21。


# NTSminer 0.1.21

## 诊断碎片离子接入默认注释流程

- **`annotate_by_spectral_match()` 新增 `validate_diagnostic_fragments = TRUE`（默认开）与
  `diagnostic_fragments_csv = NULL`**：打分结束后，自动用**包内 PFAS 专家碎片表**
  （`default_diagnostic_fragments_csv()`；Koelmel et al. 2022 + Barzen-Hanson et al. 2017，
  956 行、其中 neg 868 行）对每张实验 MS2（`@ms2data`）做诊断碎片匹配，把命中的碎片式写入
  `@ms2_annotations$diagnostic_fragments`，并将命中特征的 `@ms1_annotations$annotation_level`
  标记为 `[MS2 Verified]`。
  - 复用既有 `validate_with_ms2()`（单一实现），`batch_annotate()` / `hierarchical_annotate()`
    经 `...` 透传，随主注释流程默认生效。
  - 传 `validate_diagnostic_fragments = FALSE` 或一个不可解析的 `diagnostic_fragments_csv`
    即可关闭该步。
  - `validate_with_ms2()`（legacy 分支）在 `ms1_annotations` 缺失时由 `stop()` 改为提示并跳过
    `annotation_level` 更新（MS2 层 `diagnostic_fragments` 仍写入），使默认钩子在裸工程上也安全。
  - **注意（回归）**：默认开启会改变下游输出（`ms2_annotations` 多一列、部分 `annotation_level`
    追加 `[MS2 Verified]`），旧四表回归基线需重建。

## 置信度分级：新增 L2d（prediction − RT）

- **预测库仅 MS2 命中、但无 RT 一致证据时，此前无档可落**：`L2c`（prediction + RT）只在 RT 一致
  时给出，RT 不一致的预测命中会掉到 L4/L5。现 12 级方案升级为 **13 级**，在 `L2c` 之后新增
  **`L2d`（Probable − Prediction, no RT）**，`L3a` 及以下优先级整体顺延 +1。
  - 与 in-house（L1/L2）、online（L2a/L2b）的 ±RT 双档对齐；CFM-ID 等预测谱库无 RT 条目时也能定级。
  - `assign_confidence_levels()` 的 `prediction` 分支现为 `if (rt_any) L2c else L2d`。

## 测试

- 新增 `tests/testthat/test-diagnostic-fragments-default.R`（默认钩子命中/无 MS2 空转/可关闭）
  与 `test-confidence-levels.R` 的 L2d 用例；全量 testthat 通过。


# NTSminer 0.1.20

## 识别准确性修缮（对 `docs/analysis/NTSminer_识别准确性评估_0.1.19.md` 的定向落地）

### 置信度分级：数据库来源贯通（A1）

- **in-house 库被误判为 online，导致 L1/L2 永远不可达（P0）**：此前置信度分级靠 `Lab_ID` 的
  正则猜库类型，真实自建库（`PFAS_L02` 之类编号）既无 `db_type` 列也不匹配 `inhouse` 模式，
  于是被静默当成 online——**L1（in-house + RT）与 L2（in-house − RT）在任何数据上都不可得**。
  现打通三处：
  - `build_nts_database(..., db_type =)` 新增参数（默认 `NA`，显式 `metadata$db_type` 列优先）；
    `build_spectraverse_db()` / `build_mona_db()` 自动盖 `"online"`，`NTS_db_manager` 自建库盖 `"inhouse"`。
  - `db_type` 经 `@metadata` → 前体候选 → `@db_annotations` → 导出的主表/网络表一路带出（此前
    `annotate_by_spectral_match()` 末尾的**显式列投影**会把它丢掉）。
  - 新增导出 `set_db_type(db, db_type)` 给旧库补戳；旧库文件 `NTS_InHouse_PFAS_202511ZJ.rds`
    已就地补戳（原文件备份为 `.rds.bak`）。
  - `.nts_confidence_db_type()` 改为「显式列优先 → `Library_ID`/`ref_id` 模式兜底 → 两者皆无时
    发一次性 warning 并退回 online」，不再无声降级。

### 置信度分级：`has_rt` 语义收紧

- **`has_rt` 从「库里有 RT」改成「实测 RT 与库 RT 一致」**：原实现只要 `lib_rt` 非空即算「有 RT」，
  于是任何带 RT 的库条目都能拿到需要 RT 证据的最高档。现要求
  `rt_error_s <= rt_tol_sec`（`assign_confidence_levels(rt_tol_sec = 30)`），缺少实测 RT 误差的
  行一律不算 RT 证据。

### MS2 打分：隔离窗与隔离纯度（C4）

- **`s_iso` / `s_purity` 由恒为 `NA` 改为真实计算**：`s_iso` 由扫描的 isolation window 与实测
  前体 m/z 的偏离度得出；`s_purity` 新增 `.nts_iso_purity_for_file()`，用该文件自身的 MS1 扫描
  计算前体同位素簇纯度（`总和(前体窗内强度) / 全部强度`）。
  - 仅在 `extract_ms2(scoring = "soft"|"strict")` 时计算，默认 `scoring = "none"` 不会多读一次盘、
    基线与既有结果不变。

### MS2 质控列（C5）

- `compute_ms2_qc()` 的 QC 列（`spectral_entropy` / `top_peak_fraction` /
  `inter_scan_similarity` / `match_score_max` / `match_score_median` / `mean_purity` / `qc_grade`）
  本就随主表导出，本次只补齐 roxygen 说明：**QC 列仅作描述，不参与打分与分级决策**。

### RT 预测：分组划分 + 适用域（B5）

- **随机划分泄漏结构类似物**：`build_rt_model()` 此前按实测 RT 分层随机切分，同系物/同分异构体
  被拆到两侧，测试集误差被系统性低估。现改为**分组划分**：默认按分子指纹平均连接聚类
  （`cluster_similarity = 0.8`）成组，整组进训练或测试；也可用 `group =` 传入自定义分组
  （如化学类/骨架）。模型对象记录 `split_strategy` / `n_groups`。
- **新增适用域（applicability domain）**：建模时保存训练集指纹，`predict_rt()` 现额外返回
  `max_train_similarity`（与最近训练化合物的 Tanimoto 相似度）与 `applicability`
  （`"in_domain"` / `"out_of_domain"`，阈值 `app_min_similarity = 0.6`），把「远离训练化学空间
  却给出自信预测」的情形标出来。

### 测试与回归

- 新增 `test-db-type.R` / `test-ms2-scoring.R` / `test-rt-split.R`。
- 全量 testthat：794 通过 / 0 失败 / 0 错误（另有 2 例按设计跳过）。顺手修正两条陈旧断言：
  `test-p2p3-fixes-019.R` 的 5-4（旧 `has_rt` 语义）、`test-param-optimization.R` 的 `acq` 报错
  信息匹配（`match.arg` 抛的是 `'arg' should be one of ...`，不含参数名）。
- 用原始 PFAS MS1+MS2 数据复跑 phase 3：除新增 `db_type` 一列外，`Master_Feature_Table` /
  `Network_Nodes` / `Network_Edges` / `annotations.csv` 的所有共同列**逐值相同**（回归无副作用）；
  `db_type` 全部为 `"inhouse"`。真实库跑 `assign_confidence_levels()` 现得
  L1=55 / L2=11 / L3c=162 / L4=831 / L5=2523（此前 L1、L2 均不可达）。

# NTSminer 0.1.19

## P2/P3 修复（同一份评估报告 §14 的剩余条目，一次性清完；随 0.1.19 发布）

### P2 —— 静默失真 / 一致性

- **同名列 `Frequency` 两种含义（9-2）**：`features_long` 里的 `Frequency` 是**相对强度**
  （`intensity / int_max`，0–1），而 `mine_features()` 产出的却是**不同 NCE 的个数**（整数 ≥1）。
  现将后者改名为 `NCE_Count`，并在 `mine_features()` / `optimize_reference_database()` 的 roxygen
  里显式说明两处语义差异；相对强度列名保持不变（下游目录/测试依赖它）。
- **参考库元数据 `RT` 被硬写成 1（9-4）**：`build_spectraverse_db()` 里 `clean_meta[, RT := 1]` 会把
  "缺失 RT" 伪装成一个合法 RT=1 秒，下游 `is.na()` 判断全部失效。改为 `RT := NA_real_`。
- **`has_rt` 对字符型 `lib_rt` 判定脆弱（5-4）**：原先 `is.numeric()` 作用于整列、只返回一个逻辑值，
  于是 MSP/CSV 导入后常见的字符型 `lib_rt`（如 `"120.5"`）会让**每一行**都被判成"无 RT"。
  现改为逐元素 `as.numeric()`，数字样式的字符串仍算作 RT 证据。
- **`min_score` 默认值三方不一致（5-3）**：文档说 0.3、签名无默认、内部走 `dots$min_score` 兜底 0.3，
  而 README/vignette 例子按 0.5 理解。现统一为签名默认 `min_score = 0.5`，并显式转发到两个内部分支。
- **诊断采样扰动全局 RNG（3-2）**：`extract_spectra_from_scans()` 里用 `sample(non_zero_pmz, 3)`
  只为打印几条样例，却会推进全局随机数流、破坏用户 `set.seed()` 的可复现性。改为 `head(..., 3)`。
- **循环内 `on.exit` 句柄泄漏（3-3）**：逐文件 `on.exit(try(mzR::close(ms)), add = TRUE)` 会让句柄
  累积到函数退出才关。改为**每轮显式关闭** + 函数级兜底 `on.exit`，且关闭后置 `NULL` 防重复关。
- **局部网络扩展只走单向边（6-1）**：`build_homologue_network()` 默认 `bidirectional = FALSE` 时，
  种子引导的 BFS 只沿**存储方向**扩展，链中段的种子够不到"上游"祖先。默认改为 `TRUE`。
- **`extract_ms2` 重复读盘（4-1）**：新增可选 `xdata`，复用已读入的 `MSnbase` 实验对象
  （`filterFile` + `filterMsLevel`），读失败自动回退旧的 `readMSData(..., mode="onDisk")`，快路径不会改变结果。
- **无 precursorMz 的谱被静默丢弃（4-2）**：DIA / 缺母离子的 MS2 扫描在过滤时被无声删掉。现计数并在
  抽谱结束后 `message()` 汇报丢弃条数与原因。
- **对象校验不统一（0-1/0-2/0-3）**：`NTSProject` 新增 `validity` 与 `prototype`（裸 `new("NTSProject")`
  现在也合法），所有公开 setter 一律 `validObject()`；补上 GC 访问器缺失的 `.nts_check`，
  `get_ms2_score_grid()` 对非 S4 / 旧版缺槽对象也不再报错。
- **硬编码本机 Windows 路径（E-2）**：新增 `.nts_resolve_executable()`（env var → `PATH` → 候选路径），
  替换 `msconvert` / `PCDLManager` / Anaconda / Python 三处写死的路径；`convert_mzml()` 的 `exe` 默认改为
  `NULL` 自动探测，并在路径含空格时 `shQuote` 后再交给 `system2`。

### P3 —— 文档 / 卫生

- **清理死代码（E-4）**：删除已废弃且不可达的 `.nts_msp_*` 解析家族（16 个）、`.nts_mgf_parse_block_lines`、
  `.nts_mix_standards_cluster_ms2`、`.nts_build_eic_traces`、`.nts_is_unknown_node` / `.nts_is_anchor_node`、
  `.nts_rcdk_fp_one_smiles` / `.nts_fp_tanimoto_pair`、`.nts_ms1_from_transform_out`、`.nts_node_network_group`、
  `.nts_toxpi_join_keys`、`.nts_clean_ms2_list_to_long_dt`、`.nts_sync_ms1_from_feature_table` 等（共 24 个）。
- **删除 `R/01_db_parsers.R.bak` 构建残留（E-5）**。
- **过期示例路径（E-6）**：vignette 里写死的 `E:/DA/PFAS_screening/...` 改为占位符 `path/to/your/project/`。
- **roxygen 与实现不一致（E-7/6-2）**：`compute_fingerprints()` 更正为 base = 5 种 / **4119** 位、
  full = 8 种 / **10310** 位（原写 2378 / 8262）；`score_*` 剪枝说明更正为"移除**全部** `score_*`（含
  `score_cosine_ms2`）"，并注明节点表不带任何打分列。
- **死分支（7-4）**：`print.nts_rt_model()` 删除永不触发的 `train_metrics` 打印分支。
- **`convert_mzml` 结果并入 `@log`（E-8）**：新增可选 `project=`，传入时写一条 `convert_mzml` 日志
  （输入/输出目录、模式、解析到的 `exe`、过滤器、文件清单）并返回更新后的工程；不传时行为不变（返回文件向量）。
- **可视化归拢（8-3）**：属低影响组织性建议，评估后**不改动文件布局**（重构成本高于收益），仅在此记录结论。

新增回归测试 `tests/testthat/test-p2p3-fixes-019.R`（35 项）。全量 testthat：732 passed / 0 failed /
1 pre-existing error（`refine_xcms_params_bayes` 需 ParBayesianOptimization）/ 18 skipped。

# NTSminer 0.1.18

## P1 修复（同一份评估报告，八条一次性修完；随 0.1.18 发布）

- **`predict_rt` 用过期索引对已剪枝矩阵再剪一次（7-1）**：模型已按 `kept_colnames` 对齐并剪枝，随后又用
  `model$nzv_idx` / `cor_idx`（它们是**原始**未剪枝矩阵的列号）对 X 再剪一次，列数从 6 掉到 3，与训练宽度
  对不上（真 xgboost 直接报错）。现已删除二次剪枝，`kept_colnames` 为唯一真源。
- **`predict_logIE` 用测试集做 xgboost 早停（7-2，信息泄漏）**：早停的 `watchlist` 里放着**留出测试集**，
  等价于用测试集选轮数，测试指标因此系统性偏乐观。改为在**训练集内切 20% 内层验证集**早停，再按
  `best_iteration` 用同一批训练行重训一次；测试集只在最后评估一次。
- **`predict_logIE` 缺失描述符静默补 0（7-3）**：缺失的特征列被静默填 0 —— 而指纹位/描述符里的 0 是**真实取值**，
  于是会给出"看起来很自信"的错误预测。现在模型保存训练集列均值（`feature_centers`），缺失列按均值填补并
  **明确告警**；只有在该字段不存在的旧模型上才回退 0，且告警里说明原因。
- **丰度列识别误吞区间列（8-1）**：`.nts_abundance_colnames` 等 5 处各自维护的"保留字"列表都漏了
  `mzmin`/`mzmax`/`rtmin`/`rtmax`/`npeaks`（它们是数值型，被"取数值列"的启发式当成了样本丰度）。
  现收敛为一个共享常量 `.nts_feature_shape_cols`，五处不再各写一份、不会再漂移。
- **GNPS 导出列名不规范（8-2）**：定量表改为 FBMN 要求的字面列名 `row ID` / `row m/z` / `row retention time`，
  RT 由秒换算为**分**；MGF 由 `TITLE=` 改为 `FEATURE_ID=`（值等于定量表 `row ID`），否则 FBMN 无法把
  谱图与定量表对上。
- **参考库基线归一 / top-k 未按 NCE 分组（9-1）**：归一与 top-k 只按 `Lab_ID` 分组，导致同一化合物在低 NCE
  下采到的强峰会把高 NCE 的峰整体压到噪声阈值之下。改为 `by = c("Lab_ID", "Polarity", "NCE")`。
- **自建库极性被静默改成 NEG（3-1）**：`extract_mix_standards` 不产出 `Polarity` 列，`build_nts_database`
  于是把**整个正离子自建库**标成 NEG，正模式前筛因此拿不到任何参考且不报错。现在 `extract_mix_standards`
  落 `Polarity` 列，`build_nts_database` 用新助手 `.nts_normalize_polarity()` 容错映射（`POS`/`POSITIVE`/`+`
  → `POS` 等），**未知的非空值报错**、整列缺失时**告警**而非静默降级。
- **加合物常数漏扣电子质量（9-3）**：`nh4_m` 用的是中性 NH₄ 质量而非 NH₄⁺，`na_m`/`k_m`/`cl_m` 同样如此。
  现 `[M+Na]⁺`/`[M+K]⁺`/`[M+NH₄]⁺` 扣除电子质量、`[M+Cl]⁻` 加上电子质量（各自约 0.55 mDa）。注意 `proton`
  **不能**再扣一次：`1.00727646` 本身就是质子质量（= m(¹H) − m(e)），`[M+H]⁺`/`[M−H]⁻` 保持原语义。

`build_nts_database` 的 `metadata`、`export_to_gnps`、`extract_mix_standards` 的 roxygen 文本有更新，
`man/` 已重生成（仅这 3 个 `.Rd` 变化，`NAMESPACE` 未动）。

验证：`tests/rawdata_pfas/verify_p1_fixes_0.1.17.R`（52/52 断言通过）；
新增回归测试 `tests/testthat/test-p1-fixes-017.R`（8 个 `test_that`、51 条断言全绿）。
全量 testthat（261 个文件）：无 xgboost 时 697 passed / 0 failed / 1 error（既有）/ 18 skipped。
> 本轮修复随 **0.1.18** 发布；评估报告与验证脚本的文件名沿用当初的 0.1.17 审查标签。

# NTSminer 0.1.17

## 主流程函数修复（全流程评估遗留下来的三条高危缺陷）〔0.1.17〕

依据 `docs/analysis/NTSminer_全流程函数评估_0.1.17.md`（按数据分析流程逐阶段审查 + 实证复现）。

- **置信等级误用综合分当余弦（P0，静默给出错误分级）**：`assign_confidence_levels` 判定在线命中的
  L2a/L2b 时读取 `db$score_total`（= `0.25*mz + 0.25*rt + 0.5*cosine`），而阈值 `ms2_cosine_threshold`
  本意作用于**余弦**。因 `score_total >= score_cosine_ms2` 恒成立，该错误只会**上调**等级、永不漏判。
  现改为读取 `score_cosine_ms2`（缺失时回退 `score_cosine`；两者都没有则不构成 L2a/L2b 证据）。
  `tests/testthat/test-confidence-levels.R` 相应改用余弦列，并新增一条回归测试。
- **GC-RT 库整个写错（P0，静默产生错误库）**：`build_gc_rt_library` 里 `for (gid in unique(scores$gid)) scores[gid == gid]`
  的循环变量被 data.table 同名列遮蔽，条件恒真，**每一轮都取到全表最优参考谱**而非该化合物自己的最优谱。
  循环变量改名为 `grp`。
- **错误处理器引用未定义的 `warn_tar`（触发即中断整轮搜索）**：`search_asymmetric_smiles` 的指纹
  循环内 `tryCatch` 处理器使用 `warn_tar$n`，但 `warn_tar` 原先定义在**循环之后**，于是处理器一旦
  触发（CDK/Java 异常）便抛 `object 'warn_tar' not found`，**整批结果作废**。已把 `warn_tar <- new.env(...)`
  前移到函数层级、循环之前。（注：不可解析的 SMILES 只返回 `NULL` 并告警，不触发该处理器——原评估
  对该条触发条件的描述已更正，严重度由 P0 下调为 P1。）

三条改动均为源码内局部修改，**未改任何导出函数签名或 roxygen**，`man/` 无需重生成。
验证脚本：`outputs/verify_p0_fixes.R`（12/12 断言通过）；相关 testthat 全绿。

## 参数优化函数修复（`suggest_xcms_params` / `refine_xcms_params_bayes`）〔0.1.17〕

依据 `docs/analysis/参数优化函数评估_suggest_xcms_params_refine_xcms_params_bayes.md`（19 项 P0/P1/P2），
对特征提取前的两个参数优化函数做了系统修复。

### P0（会静默产生错误结果）

- **`obiwarp_binSize` 单位错误**：`ObiwarpParam@binSize` 是 **m/z 维度（Da）**，旧文档/示例按「秒」给出 `0.3`。
  Sonar 现返回安全默认 `0.6`（Da），并在 roxygen、`XcmsParam`、出厂模板中统一按 m/z 说明。
- **BO 目标函数可被刷分**：旧 `width/ar` 子项是自参照均匀度，对全局尺度完全不变（峰宽 ×10、强度 ×1000 分数一字不差）。
  重写 `shape` 项：对称性（0.5）+ 与理想高斯面积比的偏差（0.3）+ 峰宽占 RT 跨度的比例（0.2），**对峰宽不再不变**；
  同位素项去掉把大分子全拉到常数的 `min(0.1, …)`，改用**观测比值落在合理区间**打分。
- **BO 返回值缺少 `params`**：现同时提供顶层 CentWave 字段与 `params` 子列表，向后兼容 `refined$ppm` 与 `refined$params$ppm` 两种写法。

### P1（必然/常见崩溃或系统性科学偏差）

- **`suggest_xcms_params(sample_n = NULL)` 崩溃**（`length(files) > NULL` → `logical(0)`）已修，并新增 `seed` 参数使文件抽样可复现。
- **空 `pw_clean` 崩溃**（`if (NA < 3)`）已修，全部返回字段统一「先判 `is.finite` 再钳制」。
- **`ppm` 双重保守 / 被噪声段主导**：改用 **Tukey IQR 上界**先剔除重尾污染、再取高分位。
  真实 PFAS 标准品原先 95 分位 ≈150 ppm（恒被钳到上限 80），修后回到 ≈11 ppm。
- **`refine_xcms_params_bayes` 缺 DIA/DDA 分支**：新增 `mode = "auto"/"dda"/"dia"`，`"dia"` 只取 `function=1` 低能帧。
- **`acq` 语义澄清**：明确其是**贝叶斯采集函数**（`ei`/`ucb`/`poi`）而非 MS 采集模式；`eic_metrics` 参数误用 `isTRUE()` 的调用已修。
- **全零目标的静默退化**：若所有迭代都命中 `min_peaks` 下限，末尾发出 warning。

### P2（边界 / 性能 / 文档一致性）

- `.sonar_estimate_noise` 无拐点时回落 **中位数**（旧为 `max`，会把预过滤清空到 0% 保留）。
- 新增 `roi_gap_da` 参数（默认 0.02 Da），使 ROI 单链聚类间隔可配置并文档化其合并风险。
- 每个 ROI 的 EIC 改按 `scan_range` 相对索引，去掉 O(n_roi × n_scan) 的全长大向量分配。
- `XcmsParam` 默认值与 roxygen 对齐：`binSize` 0.015→0.008、`bw` 10→20、`min_fraction` 0.25→0.001（与 Sonar / 出厂模板一致）。
- 出厂 skill 模板补充 `obiwarp_binSize`，并接上 Sonar 的 `binSize/bw/min_fraction` 建议。

### 测试与验证

- 新增 `tests/testthat/test-param-optimization.R`（22 项）：噪声估计不变量、空分位池、ppm 重尾、BO 目标对峰宽的非不变性、同位素项高分位、参数校验。
- 抽出纯函数 `.sonar_derive_params()` / `.sonar_estimate_noise()` / `.sonar_peak_smooth()` / `.bo_peak_shape_score()` / `.bo_isotope_score()` 以便脱离 mzML 单测。
- 真实 PFAS 数据全流程复跑（`tests/rawdata_pfas/run_202511ZJ_v017/`）输出与 0.1.16 基线**逐字节一致**（仅差一处此前的非 ASCII 标签字形）——参数优化函数不进入提取/注释主链路，符合预期。

# NTSminer 0.1.16

## 全包缺陷排查修复（描述符崩溃 / `predict_rie` / `predict_logIE`）

对包做了一次系统排查（全量 `testthat` + `R CMD check`），修复以下缺陷：

### 描述符计算导致 R 进程崩溃（rcdk/CDK）

- `R/08_rt_prediction.R` 的 `.nts_compute_descriptors()` 默认取 `rcdk::get.desc.names()` 全部 53 个描述符，其中 `LongestAliphaticChainDescriptor` 在**芳环带烷基侧链**的分子（甲苯、苯甲酸、阿司匹林…）上于原生代码中**无限递归**，触发**不可捕获**的 `segfault from C stack overflow`，直接杀死 R 进程（此前 `test-logIE-*.R` 整文件崩溃即由此；`tryCatch()` 无法拦截）。
- 现在新增 `.NTS_UNSAFE_DESCRIPTORS` 黑名单并从默认集合中剔除，描述符数 53 → 52；对 RT 模型影响可忽略。

### `predict_rie()` 自动加载模型路径错误

- 旧代码用 `data/logIE_model_inhouse_{pos,neg}.rds`（**相对路径**，且漏了 `logIE_models/` 子目录），实际模型在 `data/logIE_models/`，于是**自动加载永不命中**，静默回落 LogP 近似（`docs/logIE/logIE_model_data_inventory.md` 第 6.1 条已记录）。
- 现在新增 `.nts_resolve_logIE_model_path()`，依次尝试 `system.file()`（安装包）与工作目录（源码树 / `pkgload::load_all()`），命中即加载。

### `predict_logIE()` 静默清零 POS 模型的 RT 特征

- 旧代码把 RT 元数据列**硬编码**命名为 `RT_avg`，而 POS 在库模型 `feature_cols` 用的是 `RT` → 对齐阶段 `RT` 被判缺失并**补 0**，POS 的 RT 特征恒为 0（inventory 第 6.3 条）。
- 现在由 `.nts_logIE_meta_matrix()` 按模型 `feature_cols` **决定列名**（`RT` / `RT_avg` / `MW` / `Organic_Percentage`）；不再出现「重命名→判缺→补 0」。

### 依赖与命名空间

- `xgboost` / `ranger` 为模型后端却**未在 DESCRIPTION 声明**：已加入 `Suggests`，并在 `build_logIE_model()` / `predict_logIE()` 增加 `requireNamespace()` 守护，缺失时报错清晰（此前是底层 "there is no package called 'xgboost'" 并被 `predict_rie()` 静默吞掉）。
- 移除 `ggplot2` 在 `Imports` 与 `Suggests` 的**重复声明**。
- 预测改用 S3 泛型 `predict()`，去除对 `xgboost:::predict.xgb.Booster()`（以 `:::` 触及内部符号）的依赖。

### 测试

- 新增 `tests/testthat/test-rt-descriptors.R`：黑名单存在、芳环+侧链分子不再崩溃、显式 `desc_names` 生效。
- `tests/testthat/test-logIE-concentration.R` 增补：`.nts_resolve_logIE_model_path()` 命中、`.nts_logIE_meta_matrix()` 列名遵循 `feature_cols`；对依赖 xgboost 的用例加 `skip_if_not_installed("xgboost")`。
- 全量 `testthat` **603 通过 / 0 失败 / 18 跳过**（跳过来自本机未安装的 Suggests：`xgboost`、`ranger`，以及网络 / `toxpiR` 相关用例；此前会在 `logIE` 处整体崩溃，掩盖了下游失败）。

### 打包与 `R CMD check` 清理（WARNING 8 → 0 / NOTE 5 → 2）

对上一轮 `R CMD check` 的遗留项做了一轮清理：

- **`.Rbuildignore`**：忽略 `docs/`、`.codegraph`、`plans/`、`scratch/`、`skill/`、`scripts/`、`gctof_data/`、`ai_tasks/` 及 `tests/` 下全部脚手架目录（`dist_simulate`、`rawdata_pfas`、`skill_tmp`、`feedback_fix` …）。源码包 **727 MB → 746 KB**。
- **vignette** `ntsminer-workflow.Rmd`：17 个依赖 `xlsx` 等外部包的代码块改为 `eval = FALSE`（保留 setup / 库 / 演示块）——`R CMD build` 不再硬失败，`inst/doc` 正常生成。
- **非 ASCII**：`R/*.R` 与 `man/*.Rd` 中的中文注释、全角符号、箭头、上下标统一转 ASCII（当前 **0** 处非 ASCII）。
- **依赖声明**：`bslib` / `DT` 由 `Imports` 移入 `Suggests`（未被引用）；`BiocGenerics` 补入 `Imports`；`caret` / `DBI` / `RSQLite` / `R.utils` / `signal` / `magrittr` / `xlsx` 补入 `Suggests`（消除 `::` 未声明依赖）。
- **文档补全**：为 8 个函数补齐缺失的 `@param`（`build_nts_database`、`compile_reference_database`、`build_structural_library`、`export_to_sirius_ms`、`infer_by_structure`、`precompile_target_db`、`predict_rt`、`run_sirius`）；为 `NTSProject` / `databaseClass` 补 `@name` / `@aliases`（修复 Rd 交叉引用缺失）。
- **`R/00_globals.R`（新增）**：`utils::globalVariables()` 声明 332 个 data.table NSE 列名，**消掉整段「no visible binding for global variable」NOTE**；`R/NTSminer-package.R` 补 `@importFrom stats predict median quantile sd`、`@importFrom utils tail write.table getFromNamespace`。
- **顺带修掉 4 个真实调用缺陷**：`data.table::merge()`（未导出）→ `merge()`；`stats::ordered`（不在 stats）→ `base::ordered`；`data.table::J()`（1.18 起未导出）→ `list()`；`extract_mix_standards(pol = …)` 部分参数匹配 → `polarity = …`。
- **roxygen 卫生**：`08_risk_assess.R` 中 `@noRd` 后误跟正文的块已修正；指向**无 Rd** 内部函数的 Markdown 链接（`[.nts_*]` / `[…()]` / `[invisibly]` / `[0,1]` / `[M-H]-` 等）改为代码样式，消掉「Missing link」WARNING。
- **结果**：`R CMD check`（`--no-examples --no-vignettes --no-build-vignettes --no-manual`）由 **0 ERROR / 8 WARNING / 5 NOTE → 0 ERROR / 0 WARNING / 2 NOTE**。剩余 2 NOTE 均为**本机未装 Suggests**（`xgboost`、`toxpiR`、`erah` …）所致，装齐后即消失。

## 筛选函数科学性与准确性修复（`flag_isotopes` / `flag_isf` / `recognize_isf` / `calculate_kmd`）

针对 `docs/analysis/筛选函数评估_flag_isotopes_flag_isf_calculate_kmd.md` 中列出的 P0/P1/P2 问题逐条修复。

### P0-1 `flag_isf()` / `recognize_isf()` 对同位素峰**视而不见**

- 此前二者只读 `@feature_table`，看不到 `@ms1_annotations$is_isotope`，于是**同位素峰（与母体同分子）会被当成 in-source 碎片**，甚至反过来充当母体。
- 现在：`flag_isf()` / `recognize_isf()` 新增 `exclude_isotopes = TRUE`（默认）。对 `NTSProject`，凡 `@ms1_annotations$is_isotope` 为 `TRUE` 的特征**同时**从「碎片」和「母体」两种角色中剔除；对普通表，若表内带 `is_isotope` 列则同样生效（内部亦接受 `isotope_names`）。
- 实测（3,582 特征空白数据）：ISF 候选 **200 → 139**，剔除 **61** 条由同位素峰驱动的假 ISF；修复后 `is_isf` 与 `is_isotope` 不再互相重叠。

### P0-2 卤素比例带只按**单个**卤素标定 → 多卤素同位素被漏检

- 旧 `.ISOTOPE_HALOGEN_RATIO_TABLE` 只覆盖 1 个 ³⁷Cl（[0.10, 0.50]）等单取代情形；Rdisop 理论包络 `M+2/M+0` 为 1Cl 0.32 ✓、2Cl 0.644 ✗、3Cl 0.964 ✗、4Cl 1.285 ✗、4Br 3.891 ✗ —— 凡 `n ≥ 2` 的卤素包络都被拒。
- 现在改为**理论包络拟合**：新增 `.ISOTOPE_ABUNDANCE` / `.ISOTOPE_HALOGEN_MAX_N` / `.nts_halogen_envelope_table()` / `.nts_halogen_envelope_fit()`，由天然丰度算出 `r_theory = n × minor/major`，观测 `child/parent` 与理论值相对误差 ≤ `.ISOTOPE_ENVELOPE_REL_TOL`（0.35）即接受。**旧单卤素带原样保留**作为一个非回归的「已接受子集」，故单卤素行为不变、多卤素召回增加。
- 实测：2Cl(0.644) / 3Cl(0.964) / 2S(0.089) 现被正确接受并标记 `isotope_kind = "halogen"`。（注：4Br 包络的 `M+2` 比 `M+0` 更丰，程序按「低 m/z = 子峰」的角色约定无法表达此类包络，属另一独立问题，本次未改。）

### P1-3 冲突消解残留**陈旧**同位素明细列

- `.nts_resolve_parent_conflicts()` 在剥离某行的母体角色时，此前只清 `is_isotope` / `parent_name`，残留 `isotope_kind` / `isotope_label` / `isotope_shift_da` / `isotope_c13_n`，导致 **27 行**「`is_isotope = FALSE` 但 `isotope_kind` 非空」。
- 现在剥离时**四列一并清空**（各自按列存在性守护）。

### P1-5 ISF 允许多级链（碎片同时又是母体）

- `flag_isf()` / `recognize_isf()` 在冲突归约后新增**无环约束**：凡「既是子峰又是母峰」的特征，其**母峰角色被降级**，故不会再出现 `A→B→C` 这类链式 ISF。

### P1-6 `flag_isf()` 与 `recognize_isf()` 默认值不一致

- 二者默认值统一为共享常量 `.ISF_DEFAULT_RT_TOL = 2`、`.ISF_DEFAULT_CORR_THRESH = 0.85`、`.ISF_DEFAULT_LOSS_TOL = 0.005`（此前 `recognize_isf` 的 `rt_window` / `corr_thresh` / `loss_tol` 分别为 3 / 0.9 / 0.01）。

### P2-7/8/9 其余评审项

- **P2-8 ¹⁵N/¹⁸O**：`.ISOTOPE_SHIFT_N15` / `.ISOTOPE_SHIFT_O18` 此前已定义却从未纳入检索。现在通过 `flag_isotopes(include_minor_isotopes = TRUE)` **选择性开启**（默认关闭）：开启后 ¹⁵N（0.997035）/ ¹⁸O（2.004246）作为 `isotope_kind = "minor"` 参与检索，并按 N/O 理论包络做比例校验。默认关闭是因为 ¹⁸O 与 ¹³C+2 仅差 ~0.0025 Da（< m/z 250 于 10 ppm 不可分辨）。
- **P2-9 比例校验「失败即放行」**：`.nts_validate_isotope_ratio()` 在 `child`/`parent` 强度缺失或非有限、或 `parent ≤ 0` 时**改返回 `NA`**（而非静默 `TRUE`）；`NA` 视作「未验证但接受」，并计入 `n_ratio_skipped`，函数末尾以 `message()` 报告，未验证比例不再与已验证比例混为一谈。
- **P2-7 ¹³C 上界复核**：`.ISOTOPE_C13_RATIO_MAX` **维持 0.90**（约合 80+ 个碳）。实测最大 0.878 更可能是「母体并非真单同位素峰」的数据质量信号，而非带宽问题；放宽上界会接纳更多此类误标，故不放宽，改为在代码注释中说明依据。

### 测试

- 新增 `tests/testthat/test-flag-isf.R`（7 组）：默认值统一、CF2 中性丢失配对、同位素在**碎片**与**母体**两种角色上的剔除、无环链、`recognize_isf` 与 `flag_isf` 一致。
- 新增 `tests/testthat/test-calculate-kmd.R`（7 组）：CF2 同系列打标、`KMD = NKM - KM` 自洽、同位素行排除、`min_series_size` / `md_range` / `kmd_tol` 边界。
- 扩充 `tests/testthat/test-flag-isotopes.R`：多卤素包络（2Cl / 3Cl / 2S）、包络表单调性、¹⁵N 默认不标 / 开启后标、冲突消解无陈旧明细列。
- 全量 `testthat`（除环境敏感的 `logIE-*`）**564 通过 / 0 失败**。

# NTSminer 0.1.15

## 峰面积取整、`ms1_annotations` 命名 list、MS2 全格表（`ms2_score_grid`）

### `feature_table` 峰面积取整

- `@feature_table` 的**样品强度 / 峰面积列**在写入槽位时统一**四舍五入为整数**：`set_feature_table()`、`extract_features()`、`refresh_feature_areas()` 三条路径共用 **`.nts_round_feature_intensities()`**；`mz` / `rt` / `mzmin` / `mzmax` / `rtmin` / `rtmax` 等坐标列**仍为数值型**。
- **`NA` 原样保留**（未检出的特征不会被填 0），整数列允许 `NA`。
- 面积来源列不变（仍由 `.nts_feature_area_colnames()` 判定），仅存储类型变化。

### `@ms1_annotations` 改为**命名 list**（破坏性）

- 槽位由「单个宽 `data.frame`」改为**按处理步骤命名的 `list`**。列到步骤的映射（内部 `.nts_ms1_groups()`）共 11 组、固定顺序：`base`、`eic`、`library`、`isotopes`、`isf`、`kmd`、`adduct`、`master`、`network`、`annotation`、`extra`。每个子表**自带 `name` 键**（与 `feature_table$name` 对齐），可无损重新 join；无法归类的列进 `extra`，不会被静默丢弃。
- **只存「有数据」的步骤**：槽位里**只保留真正填了值的步骤**——某步映射到的列若缺失、或整列全为 `NA`（典型如 `adduct_type` / `neutral_mass` 由 `.nts_ms1_sanitize_and_order()` 兜底建成全 `NA`、但加合物归属从未运行），该步就**不写进槽位**，于是 `names(@ms1_annotations)` 直接反映**跑过且出结果的步骤**（这套 PFAS 数据下正好 6 组：`base` / `eic` / `library` / `isotopes` / `isf` / `kmd`）。被省去的步骤**不丢信息**：其规范列在下次写入时由 `.nts_ms1_sanitize_and_order()` 重建，`get_ms1_annotations(x)` 的扁平视图**列集合与取值不变**（实测 37 列逐列 `identical`）；读取被省去的步骤（如 `get_ms1_annotations(x, "adduct")`）返回 0 行子表而非报错。裁剪在**写入时**发生，因此改造前存下的 `.rds` 仍是旧形状，下次任何一步写回（或手动 `set_ms1_annotations(x, get_ms1_annotations(x))`）即收敛。
- **`get_ms1_annotations(x)`** 仍返回**合并后的扁平 `data.frame`**（向后兼容）；**`get_ms1_annotations(x, "kmd")`** 读取单个步骤的子表；`set_ms1_annotations(x, value)` 写入扁平表（内部自动分区），**`set_ms1_annotations(x, value, "kmd")`** 只替换该步骤的列。
- 兼容**旧版保存的工程**：槽位里仍是扁平 `data.frame` 时照常读取，并在**首次写入**时升级为命名 list。
- 内部 **`.nts_ms1_sanitize_and_order()`** 契约不变（仍返回扁平 `data.frame`），约 20 处调用点无需改动；列到步骤的映射以 **`.nts_ms1_groups()`** 为唯一来源。

### 新增 `@ms2_score_grid`：MS2 匹配**全格表**

- `annotate_by_spectral_match()` 新增 **`return_score_grid = TRUE`**（默认）：把 **(库 CE) × (实验谱)** 的**每一个**打分单元格写入新槽 `@ms2_score_grid`，**包括未通过 `min_score` 的单元格**，从而可在不重新打分的前提下**离线重调 `min_score`**；`@db_annotations` 仍是**每特征 × 每化合物一行**的合并去重结果（top-1 语义不变）。
- 列：`name`、`ref_id`、`lib_CE`、`exp_spectrum`、`score_cosine`（原始点积余弦，即 `min_score` 所把关的值）、`n_matched_peaks`、`passes_min_score`、`score_ms2`（多维复合分，未过闸为 `NA`）、`score_ndp`、`ndp_raw`、`score_nl`、`score_md`、`score_diagnostic`、`is_winner`（即 `db_annotations` 对应行背后的单元格）。
- 新增导出 **`get_ms2_score_grid()`** / **`set_ms2_score_grid()`**；`extract_features()` 在重新提峰时清空该槽。**兼容 0.1.15 之前保存的工程**：旧 `.rds` 的对象根本没有该槽值（类定义后来才加），`get_ms2_score_grid()` 遇缺失槽位时返回 **0 行 `data.frame`**（即打分前状态）而不是报错；首次 `set_ms2_score_grid()` 写入时补建该槽。
- 实验谱标签补全为 **`CE<ce>#<scan_idx>@<源文件主名>`**（`db_annotations$exp_CE` 同步变化）：`scan_idx` 只在单个源文件内唯一，多文件聚合的特征此前会出现同标签谱图、导致 `is_winner` 重复；现在标签唯一，`exp_CE` 也能直接看出命中来自哪个文件。
- **`is_winner` 与最终 `@db_annotations` 保持同步**：`batch_annotate()` 会在 `annotate_by_spectral_match()` 之后**把 `@db_annotations` 去重**（每 `ref_id` 只留下 `score_total` 最优的一行）；现在它会**同步清除**那些被淘汰的特征 × 化合物对的 `is_winner` 标记（只在**缺失**的配对上置 FALSE，单元格本身**保留**，故离线重调 `min_score` 仍能看到全部 CE × 谱 的分数）。这样跑完整条流水线后 `sum(is_winner) == nrow(db_annotations)` 恒成立。

### `recognize_isf()` 中性丢失默认值与 `flag_isf()` 统一

- 中性丢失默认集合统一为 **`HF` / `H2O` / `CO2` / `CF2` / `CH2` / `SO3`**（与 `flag_isf()` 一致），`recognize_isf()` 不再自带另一套默认值。

### 测试与验证

- 新增 `tests/testthat/test-ms1-annotations-list.R`（12 组 / 86 断言）与 `tests/testthat/test-ms2-score-grid.R`（7 组 / 44 断言）：步骤分区、**只保留有数据的步骤**（全 `NA`/缺失步骤被裁掉、填充后回归、被裁步骤仍可读为 0 行）、扁平视图无损回读、旧槽位兼容、`set(step=)` 只替换该步、整数取整、全格表契约（单元格数、未过闸单元格保留、`is_winner` 唯一、与 `db_annotations` 逐格一致）、同标签谱图不再产生重复 winner、缺槽旧工程的容错读取与写入补建、`batch_annotate` 去重后 winner 数与 `db_annotations` 行数一致。
- 真实数据验证（8 针 MS1 + 2 针 MS2 PFAS；3,582 特征 / 813 个含二级谱）：`feature_table` 6 个样品列全部为整数且 `NA` 保留；`ms1_annotations` 首次写入升级为命名 list、扁平视图与旧对象逐列无差异；`@ms2_score_grid` 共 **3,082** 个单元格（过闸 2,298 / 未过闸 784）；注释阶段 `db_annotations` 84 行对应 84 个 `is_winner`，`batch_annotate` 去重到 68 行后 winner 同步收敛为 **68**，`lib_CE` / `exp_CE` / `score_cosine_ms2` / `n_matched_peaks` 与全格表逐行一致。
- **改动后端到端复跑同一套 PFAS 数据**（`run_202511ZJ_v015_int`，与改动前 `run_202511ZJ_v015` 同参数同库对比）：提峰 3,582 特征、空白 QC 后 1,690、含二级谱 813（clean_ms2 789）**与基线逐阶段一致**；注释 84 行 / 79 特征 / 68 化合物、`batch_annotate` 去重后 68 / 66 / 68，**与基线完全相同**，`score_cosine_ms2` / `score_ms2` / `score_total` / `n_matched_peaks` 逐行 `all.equal` 相等；除样品强度列由 `numeric` 变为**整数**（`new == round(基线)`、`NA` 保留、`mz`/`rt` 仍为 double）外，`feature_table` 其余列与基线**逐位相同**，`ms1_annotations` 扁平视图 37 列**无损**。


## GNPS 统一 LC-MS/MS 谱库（7 个 MSP 源）

- 新增 `tests/gnps/build_gnps_library.R`：合并 7 个 MSP 源为统一 NTSDatabase，共 **129,615 谱图**（POS 85K / NEG 45K）。
- 源构成：MONA (56K)、MASSBANK (65K)、NIST14 (5.8K)、NIH Pharm (1.5K)、MASSBANKEU (1.5K)、Pesticides POS/NEG (0.7K)。
- 产出：`data/lc_libraries/gnps_unified_library.rds`（xz 压缩）。
- 修复 `build_nts_database()` 调用：改用 `version_str` 参数（移除不存在的 `database_name` / `description`）。

## LC-MS/MS 谱库匹配：逐 CE × 逐实验谱网格

- **`annotate_by_spectral_match()`** 新增 **`ref_spectra = c("per_ce", "merged", "all")`** 参数，**默认值由 `"merged"` 改为 `"per_ce"`**：不再只与库中 `NCE == "Merged"` 的合并谱比对，而是与库中**每个 CE**（CE10/CE20/CE40）的谱图分别匹配。
- 匹配粒度改为 **(库 CE) × (实验谱)** 网格：每个 MS1 候选与「每个库 CE × 每条实验二级谱」逐一打分（`score_spectra_multidim` 多维复合分），在**余弦 ≥ `min_score`** 的组合中取**总得分最高**者；获胜的 `lib_CE` / `exp_CE` 写入 `db_annotations` **新列**。
- `min_score` 语义：只要存在任一组合余弦达标即保留该候选，报告的是达标组合中的最高总得分组合，因此 `score_cosine_ms2` 恒 **≥ `min_score`**（不会出现"复合分高但余弦不达标"的命中）。
- 一个特征可同时对应多个库化合物：`db_annotations` 保留**全部命中化合物**；分子网络 `nodes` 仍为**每特征一行**（取总得分最高化合物）。
- 兼容：显式传 `ref_spectra = "merged"` / `"all"` 可复现旧行为（单一参考谱）。
- 新增内部函数 **`.nts_ref_matrices_per_ce()`**；**`.nts_feature_ms2_peak_list()`** 支持返回 CE 标签（`ms2data` 为 `CE10#12`，`clean_ms2` 为 `CE10/CE20`）。
- 数据库匹配用 `ms2data`，分子网络仍用 `clean_ms2`。
- 验证（`projcet01_annotated_202511ZJ.rds` + 两套 PFAS 库）：命中数不变（旧 67 / 新 67，无救援、无回退）；**40/67** 命中的获胜 CE 发生变化；改按总得分排序相对按余弦排序的余弦代价中位数仅 0.001（最大 0.214）。
- 新增测试 `tests/testthat/test-annotate-per-ce-grid.R`：逐 CE 救援、网格打分、多化合物保留 + 节点单行、`batch_annotate()` 默认继承 `"per_ce"`；`tests/testthat/test-metadata-screening-mz.R` 新增 `.nts_db_ms1_table()` 的 `target_mz` 优先用例。
- **`batch_annotate()`** 不再硬编码 `ref_spectra = "merged"`（改为继承 `"per_ce"` 默认）：此前在模板里 `annotate_by_spectral_match()` 之后紧接的 `batch_annotate()` 会**回退到合并谱**、覆盖逐 CE 结果。
- 前体 m/z 解析统一**`target_mz` 优先**：`.nts_db_ms1_table()`（legacy `databaseClass` 路径）此前只用 `mz.neg`/`mz.pos`，再回退 `mass ∓ 1.007276`，从不读 `target_mz`；现改为 `target_mz` > 极性 mz > `mz` > `mass`，与 `NTSDatabase` 路径的 `.nts_metadata_screening_mz()` 一致。两套 PFAS 库中 `target_mz ≠ mass − 1.007276` 的化合物有 **11/73**（ISF/加合物/碎片离子）；旧逻辑会漏掉 **9** 个特征（FT0073、FT0151、FT0154、FT0468、FT0713、FT0857、FT1013、FT1014、FT1030）。

## GC-QTOF 管线：跨样本对齐、高分辨 refine 与 CAS 过滤谱库

### erah 检测增强（`R/00_gc_io.R`）

- **`run_erah_detection()`**：新增**强制跨样本对齐**（erah `alignComp`），新增 `cache_dir`、`compression.coef` 参数。
- m/z 精度修复：erah 谱解析器 `as.integer` → `as.numeric`，保留小数质量精度。
- 新增 **`refine_gc_spectra()`**：用 mzML 高分辨谱替换 erah 的整数 m/z 谱。

### 谱库构建

- 新增 **`build_nist_library()`**：NIST MSP → NTSDatabase（含 `NIST#` 与 `InChIKey`）。
- 新增 **`build_gc_rt_library()`**：CAS 过滤的 NIST 子库 + QTOF RT + 高分辨谱。
- NIST23 REPLIB 提取：`nist_batch_extract_replib.py`（pynist，`NISTMS_REP_LIB=4`）→ `data/gc_libraries/nist23_replib.rds`（**27,729** 化合物）。
- VOC PCDL 谱库去重：`build_gc_library()` 自动唯一化。

### 可视化

- 新增 **`extract_gc_fragment_eic()`** / **`plot_gc_fragment_eic()`**：单碎片离子 EIC 提取与可视化。

### `score_ei_spectra()` 重构（`R/06_gc_scoring.R`）

- 新增 **`mass_tol`** 参数：匹配前质量预过滤（DeepEI 风格）。
- 新增 **`filter_query`** 参数：匹配前去噪。
- 新增 **RT 过滤**（`rt_tol` / `rt_lookup` / `lib_rt_col`）与 **CAS-bridge** 匹配模式。
- 新增 **`best_only`** 选项，默认 **`FALSE`**（异构体救援，避免 winner-take-all 误判）。
- 性能：反转 m/z 索引，避免整表复制；预过滤改由 `prefilter_threshold > 0` 触发（移除 `n_lib` 门槛）。
- `build_gc_rt_library()` 修复：`prefilter_threshold` 0 → 1。

### 验证

- **BLANK-JB 加标验证**：NIST23 CAS 过滤匹配，检出率 35%、RT 准确率 83%；`data_list.xlsx` ground truth 对比。
- **Sheet2 双轨验证**（97 化合物，BLANK-JB）：NIST23 CAS 覆盖 100%（CAS 归一化修复后）；integer + refined 双轨（`best_only=FALSE`）在 cos≥0.7 达 **90%**。
- 低基峰卤代烃适用 integer 谱（保留解卷积分离），含氧物适用 refined 谱（精确 m/z 避免 O₂⁺ 干扰）；11 个化合物经 `best_only=FALSE` 救援。
- 89 测试 PASS，0 FAIL。

## 文档与仓库维护

- 新增 `docs/gc_qtof_pipeline.md` / `.html`：GC-QTOF 完整管线文档。
- 新增 `docs/ntsminer_functional_overview.md` / `.html`：双管线功能总览（含 Mermaid 图，165 个导出）。
- 仓库重组：`data/databases` 移出主仓库，清理 `R CMD check` 产物（`.check/`、`*.Rcheck/`）。


# NTSminer 0.1.14

## mzCloud LC-MS/MS 参考库集成 (Orbitrap, ESI+)

- **mzCloud Offline Library Reference 2023** 完整集成：12,889 化合物 / 573,052 MS2 谱图 → 420,190 ESI+ 谱图入库。
- 密码发现：通过 `SQLiteConnection.Changed` 静态事件拦截获得 SQLite SEE 加密密码（GUID `504025fe-1369-490b-a0c3-95e17cf3ecd1`）。
- 谱图解码：`blobMass` / `blobIntensity` 均为 8-byte double 数组，MSP 导出 1.24 GB。
- 新增 `data/mzcloud/mzcloud_NTSDatabase.rds`（545 MB xz）：POS+NEG 双极性合并，573,052 spectra。
- 化合物类别：Endogenous Metabolites (34K), Therapeutics (22K), Drugs of Abuse (21K), Pesticides/Herbicides (10K), PFAS (298)。
- 全部 Orbitrap FT 高分辨谱图，SMILES/InChIKey 覆盖率 100%。

### import_msp() 扩展

- 新增 8 个 mzCloud 特有字段解析：`COMPOUNDCLASS` → `Compound_class`, `INSTRUMENT` → `Instrument`, `MASSANALYZER` → `Mass_analyzer`, `FRAGMENTATIONMODE` → `Fragmentation_mode`, `SYNONYM` → `Synonyms`, `SPECTRUMTYPE` → `Spectrum_type`, `NEUTRALMASS` (修正 mass), `PUBCHEM` (独立行)。
- 修复 Collision energy 正则：`Collision[-_ ]energy` → `Collision[-_ ]?energy`，支持 `COLLISIONENERGY` 无分隔符格式。
- 新增 roxygen `@section mzCloud support` 和 `@examples` 示例。

### 文档

- 新增 `docs/mzcloud_integration.md`：完整提取、结构、使用文档。
- 新增 `ai_tasks/mzcloud_coverage_report.md`：化合物类别/仪器/PFAS 覆盖分析。
- 新增 `ai_tasks/convert_mzcloud_to_nts.R`：MSP → NTSDatabase 转换参考脚本。


# NTSminer 0.1.13

## GC-QTOF 数据处理管线 (EI, Agilent 7250)

> 自研 GC-EI 管线（`detect_gc_compounds()` / `run_gc_detection()` / `deconvolve_gc_spectra()` /
> `score_ei_spectra()` / `@gc_spectra`·`@gc_deconvolved` slot / Phase 16–21 迭代）见 **0.1.12** 章节；
> 本节仅记录 0.1.13 新增内容。

### erah PCA-OSD 解卷积集成

- 新增 **`run_erah_detection()`**（`R/00_gc_io.R`）：调用 erah 的 `deconvolveComp()` 进行 GC-EI 谱图解卷积，替代自研管线。在 372 个农药标品基准上达到 **64%** cos≥0.7 识别率（自研 TIC+EIC 形状聚类管线 13%）。
- 新增 **`.nts_parse_erah_spectrum()`** 内部函数：解析 erah 输出的 `"mz,int mz,int ..."` 格式谱图。
- DESCRIPTION Suggests 新增 **erah (>= 2.0.0)**，非强制依赖，用 `requireNamespace("erah")` 守卫。
- 评估脚本：`tests/erah_full_pipeline.R` — A-E 全 5 组数据 erah 解卷积 + V2 库余弦匹配。

### 谱库导入（Agilent PCDL / CEF）

- 新增 **`import_msp(polarity="EI")`**：修复极性校验，支持 EI MSP 导入；新增 `CAS#:` / `RetentionTime:` 字段提取（CAS 优先 standalone > Comments）。
- 新增 **`import_cef()`**：解析 Agilent CEF (Compound Exchange Format) XML，产出 metadata + spectra。
- **`import_pcdl()`** 文档更新：新增 SQL Server Compact 提取路径说明与 CEF 回退方案；`PCDL.cdb` 移至外部 `D:\MassHunter\PCDL\`。
- PCDL 农药谱库：`data/pcdl/pcdl_pesticides_library.rds`（1020 compounds, QTOF 级高分辨 EI）。
- **PCDL RT 校准**：线性模型（R²=0.933）将 PCDL 20 min 方法 RT 校正到实验 40 min 方法；产出 `data/pcdl/pcdl_rt_calibration_model.rds`，`data/gc_reference/GC-5mix_supplemented.xlsx` 补齐 455/473 个 RT，识别 485 个化合物。

### 测试与文档

- GC 测试合计 **89 tests, 0 失败**（test-gc-deconv.R / test-gc-ei.R / erah 集成）。
- `docs/asb_skill_collections_analysis.md`：ASB 方法学分析 (MSHub/CoreMS/mssearchr/mzapy)。
- `docs/gc_ms_methodology_comparison.md`：NTSminer vs RaMS/erah/OrgMassSpecR 对比。


# NTSminer 0.1.12

## GC-QTOF 数据处理管线（EI）

- 新增 **`suggest_xcms_params(mode = "gc")`**：GC-QTOF 专用 xcms CentWave 参数推荐（peakwidth = c(1,8), ppm = 20, snthresh = 5）。
- 新增 **`extract_gc_spectra()`**（`R/00_gc_io.R`）：从 xcms 峰表在 apex RT 提取 GC-EI 全扫描谱图，支持多扫描合并与背景扣除。
- 新增 **`read_gc_tic()`**：快速读取 GC-QTOF mzML 文件的 TIC 色谱。
- 新增 **`import_pcdl()`** / **`build_gc_library()`**：Agilent PCDL 数据库 → MSP → NTSDatabase 转换，供 GC-EI 谱库匹配使用。
- 新增 **`score_ei_spectra()`** / **`annotate_gc_by_spectral_match()`**（`R/06_gc_scoring.R`）：基于 OrgMassSpecR 余弦相似度的 EI 谱库匹配打分，结果写入 `@db_annotations`。
- 新增 **OrgMassSpecR** 依赖（Imports）：纯 R 无编译，提供 `SpectrumSimilarity()` 余弦匹配。
- NTSProject 新增 **`@gc_spectra`** slot（`named list of data.frame(mz, intensity)`），配套 `get_gc_spectra()` / `set_gc_spectra()` 访问器。
- 测试：`tests/testthat/test-gc-ei.R`。
- 新增 **`deconvolve_gc_spectra()`**：Biller-Biemann 算法 GC-EI 共流出峰解卷积，按 m/z 通道 EIC apex 聚类分离纯组分谱。NTSProject 新增 **`@gc_deconvolved`** slot，配套 `get_gc_deconvolved()` / `set_gc_deconvolved()` 访问器。
- **`score_ei_spectra()`** 新增 **`use_deconvolved`** 参数：启用时使用解卷积纯谱做库匹配，消除共流出干扰。
- 测试：`tests/testthat/test-gc-deconv.R`。
- 新增 **`detect_gc_compounds()`**：erah 式全色谱解卷积，一步完成预处理→化合物检测→纯谱提取，替代 CentWave 在 GC 路径上的使用。返回 data.table 含 compound_id、spectrum、rt、basepeak_mz 等列。
- 新增 **`run_gc_detection()`**：`detect_gc_compounds()` 的 NTSProject 集成包装，填充 `@gc_spectra`、`@feature_table`、`@ms1_annotations`。
- **Deprecated**: `suggest_xcms_params(mode = "gc")` 和 `extract_gc_spectra()` 已标记废弃，推荐使用 `detect_gc_compounds()` 替代。LC 路径不受影响。
- 新增内部辅助函数：`.nts_baseline_correct()`（running-median 基线扣除）、`.nts_smooth_eic()`（移动平均平滑）、`.nts_try_quantile()`（安全分位数）。
- **Phase 21: 修正 RT 自建 GC-EI 精准谱库 V2** — 用用户修正 RT 值的 5 组标品数据 (A-E, 372 compounds) 靶向提取谱图，匹配 GC-5mix.xlsx 获取 Formula/SMILES/CAS，构建 NTSDatabase V2 (102,620 peaks, 中位 254 ions/spectrum)。自验证中位余弦 0.875，78-84% 匹配余弦 ≥0.7。新增 `build_inhouse_library_v2.R`、`validate_inhouse_library_v2.R`、`test_ntsminer_with_v2_library.R`。产出: `data/gc_libraries/inhouse_gc_library_v2.rds` (1.0 MB)。
- **Phase 19: Agilent Ground Truth 差距分析** — 用 Agilent 靶向检出 (87/91, 95.6%) 作为 ground truth，精确评估 NTSminer 非靶向管线。NTSminer 检出 **39%** (34/87) 的 Agilent 检出化合物，主要差距在 18-24 min 高密度共流出区域。Agilent 靶向谱库自验证: 5/55 (9.1%) 余弦 ≥0.7，鉴定准确率 100%。漏检化合物中位 SNR=1,454（强信号，非低丰度），确认共流出分离是当前最大瓶颈。新增 `build_agilent_guided_library.R`、`test_gcqtof_phase19_compare.R`、`docs/phase19_gap_report.md`。
- **Phase 18: Agilent 风格峰质量评估** — 新增 `peak_ratio_min` (默认 2)、`peak_sym_min` (默认 0.6)、`delta_apex_max` (默认 5.0s)、`delta_fwhm_max` (默认 2.0) 参数。放宽进入门槛 (min_snr=2, min_ions=3, mz_bin_width=0.1)，新增后聚类质量把关：基峰对称性、簇内 apex RT 对齐、FWHM 一致性。新增 `.nts_peak_fwhm()` / `.nts_peak_symmetry()`。A500: peak_ratio_min=2 → 55 cmp, 69.2% mol ion; peak_ratio_min=3 → 28 cmp, 60.4%, 76s。测试 +10 (82 total GC)。
- **Phase 17: EIC 形状聚类解卷积** — 新增 `mz_bin_width`（默认 0.01 Da）和 `shape_min_cor`（默认 0.85）参数。用 **Pearson 相关聚类 EIC 峰形**替代纯 RT 聚类来分组碎片离子——同一化合物的所有碎片 EIC 形状高度相关 (r ≥ 0.85)，即使完全共流出也能被分离。新增内部函数 **`.nts_eic_shape_cluster()`**。A500 基准: 0.01 Da → 分子离子匹配率 **69.2%** (Phase 16: 20.9%)，谱图纯度大幅提升；0.1 Da 回归模式 **82.4%**。测试 +4 (72 total GC)。
- **`detect_gc_compounds()` Phase 16 优化**：新增 **per-m/z-bin pre-filter**（max/median 强度比 ≥ 5），仅对呈现色谱峰形态的 m/z 通道做 EIC apex 检测，过滤噪声通道。处理时间从 600s 降至 **~45s**，检出数从 45 → **54 (+20%)**。
- **`detect_gc_compounds()` 新增 `baseline_window` 参数**：可配置 running-median 基线扣除窗口宽度 (默认 31)。
- **`detect_gc_compounds()` 默认 `apex_tol` 从 1.5s 收紧到 1.0s**：更窄的 RT 容忍度避免共流出物过度合并。
- **`score_ei_spectra()` 新增 `method` 参数**：`"cosine"` (默认, OrgMassSpecR 余弦) 和 `"ei_normal"` (严格身份匹配: 查询谱每个峰必须在库谱中有 m/z 匹配, 分数=匹配强度占比)。
- 新增诊断脚本：`ai_tasks/diagnose_gc_funnel.R` 量化检测漏斗 (环 A→B→C)。
- 新增端到端验证脚本：`tests/test_gcqtof_pipeline_v3.R`。
- 测试：`tests/testthat/test-gc-deconv.R` 新增 4 个测试 (65 total GC tests)。

## logIE 预测模型（纯 rcdk，自包含）

- 新增 **`predict_logIE()`** / **`build_logIE_model()`**：基于 rcdk 指纹 + CDK 描述符 + RT/MW 的 logIE 预测（XGBoost），不依赖外部 Python/RDKit/Mordred。
- **POS 模型**：`all_fp+desc+rt_org_mw` (4617 feats), Test R²=**0.565**, `data/logIE_models/logIE_model_inhouse_pos.rds`
- **NEG 模型**：`all_fp_desc_rt_mw` (4616 feats), Test R²=**0.644**, `data/logIE_models/logIE_model_inhouse_neg.rds`
- S3 class **`nts_logIE_model`**：`print()` / `predict()` / `evaluate()` / `load` / `save` 方法。
- 新增 **`predict_rie()`** logIE 集成：SMILES → `predict_logIE()` → `rie = 10^pred_logIE`；无模型时 fallback LogP mock（向后兼容）。
- **`run_semi_quantification()`** 透传 `logIE_model` / `polarity` / `smiles_col` 参数。
- 可选的 RDKit/Mordred 描述符通过 reticulate 补充（缺失时自动降级填 0）。
- 单元测试：`tests/testthat/test-logIE-prediction.R` (11 tests)、`tests/testthat/test-logIE-concentration.R` (8 tests)。

## RT → 有机相比例插值

- 新增 **`predict_organic_percentage()`**（`R/04_feature_filter.R`）：基于 POS 梯度实测 RT↔Org% 映射，通过 `approxfun`/`splinefun` 线性插值补全有机相比例。支持 `linear`/`constant`/`spline` 三种插值，输入校验完整。
- NEG Organic_Percentage 特征测试：不提升 logIE 预测（R² 0.655→0.588），与 RT_avg 共线，保留现有模型。报告：`tests/logIEpred/logIE-ME/NEG/org_pct_test_report.md`。

## POS logIE 目录清理

- `tests/logIEpred/logIE-ME/POS/` 旧 RDKit/Mordred 描述符文件、旧 RF 子目录、旧建模文件已归档至 `archive/`；指纹文件重命名匹配 NEG 命名规范。

## SIRIUS CANOPUS 化合物分类集成

- 新增 **run_canopus** 参数于 compute_sirius_fingerprints()：设为 TRUE 时在 SIRIUS 工具链中插入 canopus 步骤。
- 解析 canopus_formula_summary.tsv：提取 NPC pathway，按 :/; 拆分为 superclass/class/subclass 四级。
- @sirius_fingerprints 新增四列：canopus_npc / canopus_superclass / canopus_class / canopus_subclass。
- 新增导出 parse_canopus_results()：对已有 batch 输出独立重解析，无需重跑 SIRIUS。

---

# NTSminer 0.1.11

## Feature table 宽窗 EIC 面积（MassHunter 对齐，可选）

- **`XcmsParam(feature_area_method =)`**：`"xcms_default"`（默认）、`"eic_rt30"`（ppm 窗 max-per-scan + RT±30 s 梯形积分）、`"apex_into"`（每样 apex 子峰 `into`）。
- **`extract_features()`**：非默认方法时样品列用新面积；`attr(feature_table, "xcms_default_intensity")` 保留 xcms `into`。
- **`refresh_feature_areas()`**：对已有 `@xdata` 重算面积，无需重跑 CentWave；实现 **`R/03_feature_area_eic.R`**。
- **`quantify_db_annotations()`**：`batch_annotate` 后按 `@db_annotations` 子集算 MH 口径 EIC 面积，**仅返回 `data.frame`**，不写回 project；新增 **`threads` / `BPPARAM`** 按样品并行（默认串行）；SnowParam 需 **`remotes::install_local()`** 安装包，不能仅靠 `load_all()`。
- 烟测 **`tests/testthat/test-feature-area-eic.R`**。

## In-house 谱库构建与谱库注释（`target_mz` / `ISF` / 无 Merged 层）

- **`build_nts_database()`** / **`optimize_spectra()`**：新增 **`spectra_mz_bin_width`**（如 `0.005` Da）、**`add_merged_layer`**（默认 `TRUE`；in-house 可 `FALSE` 仅保留 CE 层）、**`features_mz_bin_width`**（`features_long$mz_value` 分箱与按键 `max(Frequency)` 合并）。
- **`optimize_reference_database()`**：透传上述分箱/Merged 参数；`features_mz_bin_width` 在重跑 **`mine_features()`** 后二次分箱。
- **MS1 筛选 m/z 统一**：内部 **`.nts_metadata_screening_mz()`**（`R/01_db_builder.R`）优先级为 **`target_mz` → 当前极性 `mz_neg`/`mz_pos`（或 `mz.neg`/`mz.pos`）→ `Precursor_mz` → `mz` → `mass±质子`**；与 **`extract_mix_standards()`** 及 202511ZJ 建库脚本口径一致。
- **`annotate_by_spectral_match()`** / **`batch_annotate()`**（**`NTSDatabase`**）：**`.nts_meta_precursor_table()`** 默认加合子路径改用 screening m/z，不再单独用 metadata **`ISF`** 列覆盖 `mz_neg`（以 **`target_mz`** 为准）。
- **`annotate_by_spectral_match()`** / **`hierarchical_annotate()`**：实验 MS2 恢复为仅 **`@ms2data`**（逐扫描取最佳余弦）；不再在存在 **`@clean_ms2`** 时误用净化共识谱（共识谱仍供分子网络与 **`validate_with_ms2(rules_file=)`**）。
- **`annotate_by_spectral_match()`**：`@db_annotations` 按 **（特征 × `Lab_ID`）** 各保留一条最佳命中，修复同 m/z 库条目互相挤掉（原按特征只留全局最高分一条）。
- **`batch_annotate()`**（**`databaseClass`**）：**`.nts_db_ms1_table()`** 同步为 `target_mz` 优先（legacy 表仍支持 `mz.neg`/`mz.pos`）。
- **无 `NCE == "Merged"` 的谱库**：**`.nts_spectral_ref_matrix()`**（原 Merged 专用）在注释时对 per-CE 层按 **`frag_mz` 取 max intensity** 作 MS2 余弦参考谱（向后兼容仍含 Merged 的库）。
- **`mine_features()`** 前体连接：**`.nts_metadata_join_precursor()`** join **`target_mz`** 并行级解析前体 m/z。
- **测试**：**`tests/testthat/test-metadata-screening-mz.R`**。
- **外部项目（202511ZJ）**：路径与建库参数见 **`E:/DA/programs/202511ZJ/README_inhouse_database.md`**、共享配置 **`inhouse_db_config.R`**。
- **修复**：`create_inhouse_database()` 不再剥离 `target_mz`；`build_mona_db()` 聚合去重时保留 `target_mz`。

## 内标半定量流水线（CSV → 匹配 → 归一化 → RIE → 浓度）

- 新增 **`read_internal_standards()`** / **`match_internal_standards()`**：从 CSV（name, mz, rt）在 `feature_table` 中按 ppm + RT（秒）匹配内标峰。
- **`normalize_by_IS()`** 扩展 **`method`**：`geom_mean`（默认，与旧版一致）、`nearest_rt`（PFAS 长梯度推荐）、`loess`、`rt_segment`；输出 `is_anchor`、`is_normalise_method`。
- 新增 **`run_semi_quantification()`**：一键执行匹配 → 归一化 → **`predict_rie()`**（仍为 MVP mock）→ **`calculate_concentration()`**（支持 `is_match` 中 `amount` / `mean_intensity`）。
- 实现：**`R/08_internal_standards.R`**；烟测 **`tests/testthat/test-internal-standards.R`**。

## ClassyFire 化合物分类（SMILES）

- 新增 **`classify_classyfire()`**：向 [ClassyFire](http://classyfire.wishartlab.com/) 提交 SMILES（批量、tab 分隔 ID），轮询完成后返回 **Kingdom / Superclass / Class / Subclass / Direct_parent** 等 ChemOnt 层级。
- 新增 **`get_classyfire_classification()`**：按 **InChIKey** 查询已入库分类（`/entities/{key}.json`）。
- 可选 **`canonicalize = TRUE`**（默认）先调用 **`canonicalize_smiles()`**；未分类行标记 **`classyfire_status = unclassified`**。
- 实现：**`R/08_classyfire.R`**；烟测 **`tests/test05.R`**；单元测试 **`tests/testthat/test-classyfire.R`**（mock API）。
- **POST 修复**：与 **classyfireR** 一致（`POST /queries` + 原始 JSON body，非 `encode = "json"`）；默认先 **InChIKey GET**；批量 5xx 时可 **`fallback_single`** 逐条提交。
- **InChIKey lookup 修复**：移除不可靠的 `exists(".nts_toxpi_inchikeys")` 早退；lookup 始终调用同包 **rcdk** 生成 InChIKey 并记录进度；仅 **`cf_status == "found"`** 才跳过 POST；**`fallback_single = FALSE`** 且 POST 失败时 **warning** 并返回 lookup 部分结果（不再 **stop**）。新增列 **`cf_status`**。
- **InChIKey 生成（CDK 2.8+）**：**`.nts_toxpi_inchikey_from_mol()`** 改用 **`.jcall(..., "getInchiKey")`**（旧 API `getInChIKey` 在新版 CDK 下会静默失败并返回 NA）。烟测表增加 **Atrazine** 对照行（SMILES 与 CHEBI/PubChem 一致：`CCNc1nc(Cl)nc(NC(C)C)n1` → `MXWJVTOOROXGIU-...`）；InChIKey 未命中时 lookup 会再试 **canonical SMILES** 对应键。

## ISF 候选母离子推断（碎片反向搜索）

- 新增 **`infer_isf_parents()`**：对 `flag_isf` 标记的 ISF 碎片，在参考谱库（[NTSDatabase]）中按 ppm 搜索哪些化合物的 MS2 碎片 m/z 与该 ISF 碎片匹配；按 `parent_mz` 分组汇总候选化合物，按匹配碎片数降序排列。
- 输出四列写入 `@ms1_annotations`（同时填入 ISF 碎片行与母离子行）：**`isf_candidate_ids`**、**`isf_candidate_names`**、**`isf_candidate_smiles`**、**`isf_candidate_n_fragments`**（均为 `"; "` 分隔的多值）。
- 实现：**`R/04_feature_filter.R`**；单元测试 **`tests/testthat/test-infer-isf-parents.R`**（18 断言）。

---

# NTSminer 0.1.10

## ToxPi 输入预处理

- 新增 **`prepare_toxpi()`**：合并 **`predict_episuite()`** 与 **`predict_admet_risk(..., backend = "admet_ai")`** 结果，按六类 ToxPi 切片做 batch 归一化与加权；支持 `scenario`、`tox21_agg`、`bcf_baf`、`align_by`、`query_smiles`、**toxpiR** / 桌面 CSV 输出。
- 新增 **`canonicalize_smiles()`**：预测前统一 canonical SMILES（rcdk / 可选 RDKit Python），避免 EPI 改写 SMILES 导致合并错位。
- 新增 **`analyze_episuite_biodegradation()`**：EPI BIOWIN 分类 + Ready 可降解 YES/NO（对齐 `tests/toxPi_202506065.R`）。
- **合并修复**：默认 `align_by = "auto"` 在等长表时按**行序**对齐（非 full join）；结构键合并改为 inner join；检测 EPI-only / ADMET-only 行并 `stop`。
- **Biowin**：仅 **Biowin3** 为 timeframe（1–5 离散风险档）；**Biowin4/7 等为概率** `1-p`；**`Ready_biodeg_risk`** 写入 `$data` 诊断列，**不**纳入 BioDegradation 切片均值。
- **`Name` 列**：按行合并 `compound_id` / `Lab_ID` / `name` / `episuite_Name` / `query_SMILES`（修复 `episuite_Name` 存在但全 `NA` 时 `Name` 全空）；可选 **`compound_ids`** 传入筛查表标签（与 `query_smiles` 同序）。
- **输出行序**：保留合并后行序（不再按 `Name` 字母排序），与预测输入一致。
- 测试：**`tests/testthat/test-prepare-toxpi.R`**、**`tests/testthat/helper-toxpi.R`**（共享 fixture）。

## ToxPi 打分与作图（toxpiR）

- 新增 **`run_toxpi()`**：由 **`prepare_toxpi()`** 的切片级得分（0–1）构建 [toxpiR](https://cran.r-project.org/package=toxpiR) **`TxpModel`** 并调用 **`txpCalculateScores`**；可选 **`extra_slices`**（如 **`Ready_biodeg_risk`**）、**`weight_multiplier`**（默认 100）。
- 新增 **`plot_toxpi()`**、**`toxpi_scores()`**、**`toxpi_slice_scores()`**：封装 toxpiR 绘图与得分提取；**`toxpi_scores()`** 以 **`TxpResult@txpIDs`** 命名（toxpiR 原始 **`txpScores`** 无 `names`）。
- **`DESCRIPTION` `Suggests: toxpiR`**；测试：**`tests/testthat/test-run-toxpi.R`**（`skip_if_not_installed("toxpiR")`）。
- 文档：**`tests/toxpred/README.md`**、**`AGENTS.md`** 步骤 11、**`README.md`** ToxPi 小节。
- 烟测脚本：**`tests/test04.R`**（分 Phase 0–5，输出 **`tests/toxpred/output/test04_<timestamp>/`**）。

---

# NTSminer 0.1.9

## 移除 OPERA 集成

- **移除** **`predict_opera()`** 及全部 OPERA CLI 辅助函数：OPERA 2.9 依赖 **MATLAB Runtime 9.12（~2 GB）**，无公网 API，部署成本过高。
- 环境归趋 / ADMET 仍使用 **EPI Suite**（`predict_properties` / `predict_episuite`）与 **ADMETlab 3.0 + ADMET-AI**（`predict_admet_risk`）。
- 烟测 **`tests/toxpred/`** 同步删除 OPERA 开关与文档。

---

# NTSminer 0.1.8

## OPERA QSAR 预测（本地 CLI）— 已在 0.1.9 移除

- 新增 **`predict_opera()`**：调用 NIEHS/EPA **OPERA 2.9** 命令行（SMILES → `.smi` 批量 → CSV），覆盖 logP、BCF、Biodeg、MP/BP/VP/WS、Koc、HL、Ready_biodeg、CATMoS 等 QSAR endpoint；需本地安装并设置 **`NTSminer_OPERA_EXE`**。
- 烟测：**`tests/toxpred/run_toxpred_smoke.R`** 增加 `run_opera` 开关（默认 `FALSE`）。

## ADMET-AI 解析修复（0.1.7 后续）

- **`predict_admet_ai()`** / **`predict_admet_risk()`**：默认 **reticulate**（与 **`tests/toxPi_202506065.R`** 一致）；`model$predict()` 返回 R `list` 时用 **`do.call(cbind, ...)`** 解析；支持 **`NTSminer_ADMET_AI_PYTHON`**。

---

# NTSminer 0.1.7

## 环境归趋 / ADMET 预测（EPI Suite + ADMETlab 3.0 + ADMET-AI）

- **`predict_properties()`** 现对接 **EPI Suite**（`episuite.dev`）：支持 **`NTSProject`**（从 `master_smiles` / `SMILES` 等解析）、SMILES 字符向量或含 `SMILES` 的 `data.frame`；宽表写入 **`@physchem_proxy`**，并保留 **`LogP_pred`**（来自 `LogKow_est`）等兼容列。实现：**`R/08_toxpred_engines.R`** + **`R/08_risk_assess.R`**。
- **`predict_admet_risk()`** 对接 **ADMETlab 3.0** REST 与/或本地 **ADMET-AI**（`backend = "admetlab3" | "admet_ai" | "both"`）。
- 新增导出：**`predict_episuite()`**、**`predict_admetlab3()`**、**`predict_admet_ai()`**。
- 测试：**`tests/testthat/test-toxpred-engines.R`**（httr mock + 解析 fixture）；外部烟测仍用 **`tests/toxpred/run_toxpred_smoke.R`**（现调用包内 API）。
- **`Suggests`** 增加 **`reticulate`**（ADMET-AI Python 路径）。

---

# NTSminer 0.1.6

## 多源注释合并（谱库 + 分子网络）

- 新增 **`merge_annotation_layers()`**：将 **`batch_annotate`** 的谱库命中与可选 **`annotate_by_dual_network`** 结果（**`net_frag_ann` / `net_nl_ann` / `net_md_ann` / `net_homologue_ann`**）按特征 **`feature`** 合并为宽表；**`annotation_method`**、**`ref_id`**、**`compound_name`**、**`smiles`** 等多值以 **`"; "`** 去重拼接；并生成 **`ref_id_Spectral_DB`**、**`ref_id_Network_DF`** 等分列便于查冲突。
- 新增 **`count_annotation_layers()`**：统计各层（**`Spectral_DB`**、**`Network_DF`**、**`Network_NL`**、**`Network_MD`**、**`Homologue`**）有 **`ref_id`** 的特征数及网络种子/传播层计数。
- 实现文件：**`R/07_annotation_merge.R`**；**`tests/testthat/test-merge-annotation-layers.R`**。

## MS2 分子网络：`max_frag_freq` 与结构传播

- **`max_frag_freq`**（**`generate_network_data`** / **`get_network_data`**）：修正为对 **`@clean_ms2`** 共识峰 **`frequency` 的下限（floor）**——保留 **`frequency >= max_frag_freq`**（及缺失/非有限 `frequency` 的峰），不再使用错误的 **`frequency <= max`**。默认仍为 **0.8**；筛查常用更严值（如 **0.3**）须在调用时显式传入。
- **`annotate_by_dual_network()`**：传播成功时除 **`propagated_*`** 外，同步 **`ref_id`、`Compound.name`、`SMILES`**（及节点上 **`Library_ID` / `Annotation`** 若存在）。负模式中性质量：**`target_neutral = mz + 1.007276`**（实现为 **`mz - mass_offset`**，`mass_offset = -1.007276`）。
- 新增测试：**`test-max-frag-freq-floor.R`**、**`test-annotate-dual-network-identity.R`**。

## `tests/test03.R`：XCMS 串联优化 + 三网注释（DDA）

- **Phase 0**：**`suggest_xcms_params(..., mode = "dda")`** → 可选 **`refine_xcms_params_bayes()`**（**`run_xcms_bo_refine`**，默认 **`TRUE`**；需 **`ParBayesianOptimization`**）→ **`XcmsParam`**（CentWave 来自 BO/Sonar，分组参数来自 Sonar）。
- **推荐顺序**：**`batch_annotate`** → **`generate_network_data`**（可按 **`ms2_edge_types`** 分别建 DF/NL/MD 网）→ **`annotate_by_dual_network`** → **`merge_annotation_layers`** / **`count_annotation_layers`**。
- 同系物部分仍为 **Step 1–3**（**`get_network_data(homologue)`** → **`build_homologue_network`** → **`infer_homologue_structures`**）。

## 开发者文档

- **`handoff.md` / `project_memory.md` / `AGENTS.md`**：分子网络注释分层、**`batch_annotate` 须在 `generate_network_data` 之前**、谱库 vs 网络传播差异、XCMS 优化宜用 **基质加标/QC**（非单独空白）等已记入交接与踩坑日志。

---

# NTSminer 0.1.5

## 文库子类 × 碎片特征：多分辨率聚类

- **`summarize_subclass_fragment_features()`**：新增 **`mz_resolutions`**（**命名 `list`**）。每项为 **`list(digits = …)`**（与 **`mz_digits`** 同形，可 **`c(DF=3, NL=3, MD=3)`**）或 **`list(ppm = …)`**（可选 **`da`**）。一次运行可叠加多套 m/z 容差；**`by_feature`** 与 **`subclass_coverage`** 增加 **`resolution_label`** 列。不传 **`mz_resolutions`** 时行为与旧版单容差参数（**`mz_tol_ppm` / `mz_digits`**）一致，且无 **`resolution_label`**。
- **`plot_subclass_fragment_features()`**：新增 **`resolution_label`**，从多分辨率 **`$by_feature`** 中择一绘图；若 **`resolution_label` 为 `NULL`** 且表中存在该列，默认取第一个分辨率。

```r
sumry <- summarize_subclass_fragment_features(db, class_col = "Second_Class",
  mz_resolutions = list(
    "digits3" = list(digits = c(DF = 3, NL = 3, MD = 3)),
    "digits2" = list(digits = c(DF = 2, NL = 2, MD = 2)),
    "ppm50"   = list(ppm = 50)
  ))
plot_subclass_fragment_features(sumry, resolution_label = "digits3")
```

## SIRIUS 指纹生成

- 新增 **`compute_sirius_fingerprints()`**：连接 SIRIUS 6 CLI，将 `@clean_ms2` 按 batch 导出为 `.ms` 文件，调用 `formulas + fingerprints + write-summaries` 工具链，从各 batch 的 `formula_identifications.tsv` 中提取分子式、`SiriusScore`、`numExplainedPeaks`、`formulaRank` 等，合并后存入 `@sirius_fingerprints`。
- 新增 **`@sirius_fingerprints` slot**（`data.table`）及配套 **`get_sirius_fingerprints()`** / **`set_sirius_fingerprints()`** 访问器。每行对应一个 feature 的 top 分子式候选，含 `name`、`adduct`、`molecular_formula`、`rank`、`sirius_score`、`zodiac_score`、`confidence`、`explained_peaks`、`fingerprint`（空占位）、`fp_version` 列。
- **实现细节与已知约束**：
  - 多 batch 写入独立子目录（`batch_1/`…`batch_N/`），避免 TSV 覆写。
  - 列名自动检测（大小写不敏感 + 模糊关键词匹配），排除 `SiriusScoreNormalized` 等变体。
  - `mappingFeatureId` 带 `trimws()` 去除尾部空格，并做 suffix 匹配回 `@clean_ms2` 的原始 feature 名。
  - SIRIUS 6 不通过 `write-summaries` 导出 `.fpt` 指纹文件；`fingerprint` 列保留为空 `numeric(0)`，待未来通过 API 实现。
  - `zodiac_score` 和 `confidence` 在当前公式注释摘要中不可得（`ZodiacScore` 需单独用 `zodiac` 工具链；`Confidence` 在 `structure_identifications` 摘要中）。
- **实测结果**（test02.R，227 features → 214 行）：`SiriusScore` / `explained_peaks` / `rank` 100% 非 NA，列映射与批次合并通过验证。
- **修复**：成功写入 `@sirius_fingerprints` 后，日志与汇总信息误用未定义变量 `out_df` 导致报错；已改为使用实际结果表 **`out`**。

## SIRIUS 指纹向量 API 提取

- 新增 **`extract_sirius_fingerprints_api()`**：连接 SIRIUS 6.x REST API（自动启停或接入已有服务），从 `compute_sirius_fingerprints()` 输出的 `.sirius` 项目文件中提取真实 CSI:FingerID 指纹向量，填充 `@sirius_fingerprints` 的 `fingerprint` 列。
- **实现方式**：
  - 自动定位 `.sirius` 项目目录（`sirius_projects_dir` / `output_dir` / 搜索 `tempdir()`）。
  - 收集项目至临时工作区，启动 `sirius service --port XXXX --headless`（或复用 `api_base` 指定的已有服务）。
  - 通过 REST API 扫描项目 → 特征 → top-N 公式候选（`optFields=predictedFingerprint` 或无此字段时回退专用指纹端点）。
  - 遍历 `@sirius_fingerprints` 的每行，将映射到的指纹向量写入 `fingerprint` 列的 `numeric()` 元素。
- **依赖变更**：`httr` 和 `jsonlite` 从 `Suggests` 移至 `Imports`。
- **已知约束**：需要 SIRIUS 账号已登录（CSI:FingerID 功能）；`.sirius` 文件在调用前仍存于磁盘。

## SIRIUS 函数健壮性修复（审阅跟进，0.1.5 内）

- **`compute_sirius_fingerprints()`**：
  - **`repair_id()`** 特征名修复：`endsWith` 后缀匹配改为带边界正则 `(^|_)(feature_names…)$`，排除 `FT01` 误匹配非后缀（如 `batch_FT011`）。
  - **CSV 列名检测**：改为精确匹配优先（`^SiriusScore$` / `^ZodiacScore$` / `^explainedIntensity$` / `^rank$`），模糊关键词仅作回退。
  - **必需列校验**：adduct / molecularFormula 列缺失时 `stop` 并打印可用列名，避免静默产生空列。
- **`extract_sirius_fingerprints_api()`**：
  - **`.kill_all_sirius_services()`**：移除 `taskkill /f /im java.exe`（会杀死用户全部 Java 进程）。无端口参数时，改为遍历 `LISTENING` 端口并通过 `tasklist` 过滤仅杀 `java.exe` 进程。
  - **代理旁路**：`no_proxy = "*"` → `no_proxy = "…,localhost,127.0.0.1"`，保留用户原有代理设置。
  - **端口正则**：`:%d\\s+` → `:%d\\b`，防止端口 80 误匹配 `:8080`。
  - **指纹跳过日志**：API 失败/缺指纹/无映射的特征名收集到 `failed_names`，结束时汇总输出。
  - **进度条守卫**：`on.exit(close(pb))` 改为 `if (!is.null(pb)) close(pb)`；加 `anyNA` 校验指纹向量。

## @sirius_fingerprints 数据表重构（0.1.5 内）

- **槽位类型**：`@sirius_fingerprints` 从 `data.frame` 改为 `data.table`，与 `NTSDatabase` 的 `metadata`/`spectra_long`/`features_long` 槽位类型一致。
- **setter 校验**：`set_sirius_fingerprints()` 新增列名校验（`.validate_sirius_fingerprints`），写入时检查 10 个必需列及 `fingerprint` 是否为 list-column。
- **原型初始化**：`create_nts_project()` 中 `@sirius_fingerprints` 原型从 `data.frame()` 改为 `data.table::data.table()`。
- **去除 `as.data.frame()` 转换**：`compute_sirius_fingerprints()` 末尾不再将 `data.table` 转为 `data.frame`，直接写入槽位。

## parse_sirius_fingerprints 统一写入（0.1.5 内）

- **`parse_sirius_fingerprints()`**（SIRIUS 5.x `.fpt` 路径）输出位置从 `@network_results$sirius_fingerprints`（list）改为 **`@sirius_fingerprints`**（`data.table`，10 列标准 schema）。
- **无元数据声明**：`.fpt` 文件仅含原始指纹向量，不含 adduct / molecular_formula / scores；非指纹列全部为 `NA`。roxygen `@details` 说明迁移至 `compute_sirius_fingerprints()` + `extract_sirius_fingerprints_api()` 可获取完整元数据。
- **文档同步**：`AGENTS.md` 速查新增"外挂引擎 / SIRIUS 指纹"分类，调用顺序中插入可选 SIRIUS 步骤；`man/parse_sirius_fingerprints.Rd` 与 `man/get_sirius_fingerprints.Rd` 随 roxygen 重生成。

## MoNA / MSP 导入（0.1.5 内）

- **`import_msp()`**：支持 **`PRECURSORTYPE:`** / **`PrecursorType:`**（连写）与 **`IONMODE:`** 等离子模式行；此前仅识别 `Precursor_type:` 等分隔写法，导致 `Precursor_type` 全为 `NA`，**`build_mona_db()`** 在「标准加合离子类型」步骤清空全表。
- **`build_mona_db()`**：若所有记录的 **`Instrument_type`** 均缺失，则不再套用 HRMS 仪器名正则（避免无仪器字段的整理版 MSP 被整表丢弃）。

## 同系物网络 API（0.1.5 内）

- **`build_homologue_network()`**：新增可选参数 **`network_data`**，复用已有 **`get_network_data(type = "homologue")`** 结果，避免重复构建 KMD 边；**`infer_homologue_structures()`** 仍应使用本函数输出（边上需 **`transformation`**），而非裸 KMD 表。
- **`tests/test03.R`**：演示三步管线与 `network_data` 复用。

---

# NTSminer 0.1.4

## Breaking changes

- 移除已长期 **`.Deprecated()`** 的 PFAS 名称别名（请改用通用函数）：
  - `compile_pfas_database` → `compile_reference_database`
  - `optimize_pfas_database` → `optimize_reference_database`
  - `build_pfas_structural_library` → `build_structural_library`
- `adjust_pfas_formula` / `adjust_pfas_formula_v` **仍保留**：仅为 `adjust_homologue_formula(..., base_formula = "CF2")` 的便捷封装，测试与文档仍依赖；若删除需同步改写测试与所有调用方。
- **`purify_ms2_spectra()`** 默认 **`merge_strategy = "ce_top_intensity"`**（多 CE 数据下 **`spectrum`** 与旧版 **`union_then_filter`** 可能不同）；回退请显式传 **`merge_strategy = "union_then_filter"`**（详见本节 DDA / MS2）。

## DDA / MS2 管线增强（P1–P4）

- **`extract_ms2()`**：每个 `@ms2data` scan 元素无条件写入 **`ce` / `polarity` / `isolation_lo` / `isolation_hi`**（缺失时为 `NA`）；新增 **`scoring`（`none`/`soft`/`strict`）**、**`weights`**、**`min_match_score`**、**`isolation_width_da`**。默认 **`scoring = "none"`** 与旧版行为一致；`soft`/`strict` 写入 **`match_score`** 与子分项（**`s_iso` / `s_purity`** 需后续接入 MS1，当前为 `NA`）。
- **`purify_ms2_spectra()`**：默认 **`merge_strategy = "ce_top_intensity"`** —— 当存在多档 **CE** 时，总 **`spectrum`** 按 m/z 对齐后取各 CE 子共识中的**最强强度**，碎片可见性可能优于旧版「先合并再过滤」；需要与旧版数值完全一致时使用 **`merge_strategy = "union_then_filter"`**。新增 **`by_ce`**、**`ce_bin`**、**`keep_ce_subspectra`**；`@clean_ms2[[k]]` 可含 **`ce_levels` / `ce_subspectra` / `ce_summary`**。**`@ms2_annotations`** 同步 **`n_ce_levels` / `ce_min` / `ce_max`**。
- **precursor 邻域清理**：新增 **`precursor_filter`**（命名 list）及导出 **`precursor_filter_template()`**（`default` / `qtof` / `orbitrap` / `permissive` / `strict`）。**`precursor_filter = NULL`** 时仍为 **`neutral_loss >= 1 Da`**，与旧版一致。
- **MS2 QC**：新增 **`compute_ms2_qc()`**、**`ms2_qc_thresholds()`**，向 **`@ms2_annotations`** 追加 **`spectral_entropy`**、**`top_peak_fraction`**、**`inter_scan_similarity`**、**`match_score_*`**、**`mean_purity`**、**`qc_grade`**；**不自动**改变 **`score_spectra_multidim`** 等打分权重。
- **DESCRIPTION**：**`Imports:`** 增补 **`methods`**（与 **S4 / NAMESPACE `exportClasses`** 一致，便于 **`R CMD check`** 通过依赖扫描）。

## MS1 / EIC 形状与绘图（`extract_features`、`@xdata`、`plot_feature_eic`）

- **`extract_features()` + `XcmsParam(eic_metrics = TRUE)`**：NTSProject 分支可向 **`@ms1_annotations`** 写入 4 个基于 **`xcms::chromPeaks()`** 的 EIC 形状代理列。若内部未附着 xcms 对象（`.xdata` 缺失）会 **message** 提示；**`refresh_eic_shape_metrics(project)`** 可在 **`project@xdata`** 仍可用时从 **`featureChromatograms` 管线结果**中补写或更新这 4 列，无需重跑特征提取。
- **`@ms1_annotations` 列序**：**`.nts_ms1_sanitize_and_order()`** 将上述 4 列固定在 **`mz` / `rt` 之后、同位素与 KMD 等 tail 列之前**，便于查看与导出。
- **`plot_feature_eic(project, feature_name, …)`**（依赖 **`project@xdata`**，与 **`extract_features`** 中 **`findChromPeaks → adjustRtime → groupChromPeaks → fillChromPeaks`** 后的对象一致，含 RT 对齐与 gap filling）：
  - **提速**：默认优先 **`xcms::featureChromatograms()`**（紧致 m/z 窗口 + **`expandRt = rt_tol`**，一次调用覆盖所有 trace × 样品）；否则 **`xcms::chromatogram(mz` 矩阵, `rt` 并集窗口)**。提取时强制 **`BiocParallel::SerialParam()`**，减轻 Windows 上多进程 **`SnowParam`** 在 chromatogram 上的额外开销。
  - **同位素**：默认 **`show_isotopes = TRUE`**，按 **`flag_isotopes`** 在 **`@ms1_annotations`** 中指到同一 **`parent_name`** 的子峰叠画；可选 **`include_parent_if_child`**。**`use_bounds`** 控制是否优先用 **`mzmin/mzmax`、`rtmin/rtmax`**（否则 **`mz_tol` / `rt_tol`** 对称窗）；默认 **`rt_tol = 15`**（秒）。
  - **`graphics`**：**`"ggplot"`**（默认）：多样品同面板叠加，**颜色 = 样品**、**线型 = trace**（M vs 同位素链）；可选 **`show_integration`** 用 **`rtmin/rtmax`** 画半透明 ribbon；**`smooth`**（`moving` / `savgol` / `none`）与 **`smooth_window`** 仅作用于 ggplot 路径。
  - **`graphics = "xcms"`**：调用 **xcms** 对 **`XChromatograms`** 的 **`plot()`**，按 vignette 习惯设 **`peakBg = colours[chromPeaks(object)[, "sample"]]`** 以显示 **chromPeak 积分**；**`peak_type`** 对应 **`peakType`**；**`sample_col_map`** 为按 **`basename(fileNames)`** 命名的颜色向量；**`...`** 传入 **`plot()`**。此路径需要 **`featureChromatograms`** 得到 **`XChromatograms`**（即图中每条 trace 必须都是 **`featureDefinitions(xdata)`** 里的 feature）；若仅回落到无 **`chromPeaks` 附加**的 **`MChromatograms`**，会报错并提示改用 **`graphics = "ggplot"`** 或保证同位素峰均已 **group** 为 feature。
  - **`engine`**：**`auto` / `feature` / `chromatogram`**，控制是否强制 **`featureChromatograms`** 或 ROI **`chromatogram`**。

## 代码健壮性与 API 小补（审阅跟进，仍属 0.1.4）

- **`XcmsParam()`**：在构造阶段校验 **`threads` / `ppm` / `peakwidth` / `snthresh` / `noise` / `prefilter` / `integrate` / `mzdiff` / `scanrange` / `binSize` / `bw` / `min_fraction` / `obiwarp_binSize`** 等基本类型与范围，避免错误参数拖到 **`xcms::CentWaveParam()`** 才失败。
- **`calculate_exact_mass()`**：去掉对原始 **`formula`** 向量长度的误检，仅在 **`as.character(.)[1L]`** 后判断空 / **NA**。
- **`predict_admet_risk()`**：文档明确**无远程后端、预测列恒为 **NA****；去掉无意义的 **`tryCatch`** 尾；**`timeout_sec`** 仍做规范化（供将来 HTTP 使用）。
- **`create_inhouse_database()`**：MS2 片段列表**预分配**再按索引写入，避免 **`list()` 逐条增长**。
- **`NTSProject` 访问器**：新增 **`get_db_annotations` / `set_db_annotations`**、**`get_physchem_proxy` / `set_physchem_proxy`**、**`get_ms2data` / `set_ms2data`**、**`get_clean_ms2` / `set_clean_ms2`**、**`get_project_log` / `set_project_log`**。
- **测试**：**`tests/testthat/test-review-hardening.R`**（**`XcmsParam`**、**`calculate_exact_mass`**、访问器烟测）；**`tests/test01.R` 等**仍为人工/数据驱动脚本，**paramounter** 大脚本已删除且不再恢复。
- **`extract_ms2()`**：当 **`@ms2data` 顶层 feature 分组数** 多于 **`ms1_annotations`** 行数时改为 **`warning`**（不再 **`stop`**），以便在对齐 MS1/MS2 名称时可继续排查下游。
- **`purify_ms2_spectra()`**：按 feature 的 **`clean_dt` 子集** 由循环内全表扫描改为 **`split()` 预分组**，大表场景下减少冗余扫描。
- **`run_sirius()`**：默认 **`sirius_exe = NULL`**，依次尝试环境变量 **`NTSminer_SIRIUS_EXE`**、**`SIRIUS_EXE`**；去掉硬编码本机路径。
- **Shiny**：**`inst/shiny_app/app.R`** 标题栏版本号改为 **`packageVersion("NTSminer")`**，与 **DESCRIPTION** 一致。

## 文库子类 × 碎片特征（`NTSDatabase`）

- **`summarize_subclass_fragment_features(db, class_col = "Second_Class", …)`**：按 **`Lab_ID`** 对齐 **`@metadata`** 与 **`@features_long`**，按子类 × **`Feature_Type`** 对 m/z 成簇（默认 **`mz_tol_ppm`**，或 **`mz_tol_ppm = NULL`** + **`mz_digits`**）。返回 **`by_feature`**、**`class_counts`**、**`subclass_coverage`**（含 **`max_prevalence`**）、**`params`**；可选 **`prevalence_denom = "with_features"`**、**`quiet`**。**多套容差（`mz_resolutions`）与绘图 `resolution_label` 见 0.1.5。**
- **`plot_subclass_fragment_features(summary, …)`**：子类 × 碎片普及率热图；**`subclasses`** 筛选时若子类仅因 **`min_prevalence`** 未进入 **`by_feature`**，对完整 **`summary` list** 会给出说明。

## 文档与开发者说明

- **`extract_features()` + `XcmsParam(eic_metrics = TRUE)`** 与 **EIC 形状列 / `refresh_eic_shape_metrics` / `plot_feature_eic`**：见上文 **「MS1 / EIC 形状与绘图」**；实现细节（**`peakidx`** 容错、**`.nts_apply_eic_shape_to_ms1_df()`**）见 **`R/03_feature_extract.R`**；EIC 列名常量 **`.nts_ms1_eic_shape_colnames()`** 位于 **`R/00_eic_shape_constants.R`**（**`R/02_project_core.R`** 的列序逻辑引用之）。
- **`refine_xcms_params_bayes()`**：新增**可选**的贝叶斯优化接口（弱依赖 **`ParBayesianOptimization`**，加入 **`Suggests`**）。在 1--3 个 QC 文件上对 CentWave 8 维参数（`peakwidth`、`ppm`、`mzdiff`、`snthresh`、`noise`、`prefilter` 双项）做 GP-BO 搜索；目标 = **峰形分 + 简易 M+1 同位素分** 加权（脚本源自 **`tests/20260405_01_ParBayes_Pick_optimize.R`**）。返回值的 CentWave 字段与 **`suggest_xcms_params()`** 同形，便于互换塞入 **`XcmsParam()`**；**不优化** PeakDensity / Obiwarp（沿用 `suggest_xcms_params()` 或手动设）。典型耗时 20--90 分钟，远高于 Sonar；适合 Sonar 初值不理想时作为 refine。
- **`suggest_xcms_params()`**：
  - 除 CentWave 外，增加 **`mzdiff`**（默认约束在 **[-0.001, 0.01]**，由 ROI 内 m/z 跨度推导）、**`prefilter`** 扫描数/强度上下限（默认与典型 BO 窗口一致，可调），以及对 **`PeakDensity`**（**`min_fraction` / `binSize` / `bw`**）与 **`Obiwarp`**（**`obiwarp_binSize`**）的粗建议；**`XcmsParam()`** / **`extract_features()`** 支持 **`obiwarp_binSize`**（`NULL` 时仍为 **0.6**）。
  - 新增 **`mode = c("auto","dda","dia")`**：`"dia"` 针对 **Waters MSE** mzML（整文件 `msLevel=1`，低/高能由 `spectrumId` 的 `function=N` 区分），仅取 **`function=1`** 低能 scan 参与 Sonar，避免高能碎片污染 `ppm` / `peakwidth`；`"auto"`（默认）按 spectrumId 自动判定。
  - 新增 **`min_fraction_default = 0.1`**：Sonar 只用每批 QC 文件无法推断真实平行设计，因此 `min_fraction` 改为用户可调的占位默认（NTS 无平行场景下 0.1 合理），并在末尾 `message` 中提示按实验设计在 **`XcmsParam()`** 覆写。
- **`filter_blanks()`** 的实现与 roxygen 位于 **`R/04_feature_filter.R`**（与 MS1 筛选同文件）；**`R/03_feature_extract.R`** 中导出顺序为 **`XcmsParam` → `extract_features` → `convert_mzml`**；**`export_mirror_report()`** 位于 **`R/09_plotting.R`** 末尾。参见根目录 **`AGENTS.md`**「源码分层」与「典型调用顺序」。
- **`tests/test02.R`**：在已有 centroid mzML（`data/PFAS/mzML/MS1` + `MS2`）上跑通 **`suggest_xcms_params` → MS1 → `extract_ms2`/`purify_ms2_spectra` → 可选 `validate_with_ms2` → `compute_ms2_qc`**。不依赖全库编译或本机 xlsx 路径。
- **`scratch/`**：示例脚本 **`run_pipeline_from_suggest_xcms.R`**（从参数建议到 MS2 QC）、**`test_p2_ce_consensus.R`**、**`test_p3_p4_ms2qc.R`** 等与上述管线一致，供本地烟测。

# NTSminer 0.1.2

## Major Improvements & Optimizations

* **Molecular network**: Optimized molecular network generation.
* **MCS / propagation**: Enhanced MCS propagation logic for targeted / network annotation workflows.
* **Node metadata**: Cleaned up node metadata in network outputs.

# NTSminer 0.1.1

## Major Improvements & Optimizations
*   **Performance Boost**: Refactored `build_structural_library` to use chunked batch processing for chemical fingerprint generation. Achieved a ~40x speedup for substructure fingerprints.
*   **Stability**: Fixed JNI/JVM crashes on Windows by removing incompatible timeout wrappers in `rcdk` calls.
*   **Targeted Propagation**: Introduced `annotate_neighbors_from_csv()` for high-efficiency structure propagation. This engine filters candidates by MS1 mass before computing fingerprints, making it suitable for very large external libraries.
*   **Enhanced Scoring**: Added Maximum Common Substructure (MCS) scoring as an optional verification layer for network annotations.
*   **BFS Logic**: Aligned targeted propagation with breadth-first search (BFS) principles, including `Propagation_Level` tracking and MS2 edge-weight based competition.
*   **Memory Efficiency**: Implemented hash-based caching for molecules and fingerprints to minimize redundant parsing of SMILES strings.

# NTSminer 0.1.0
*   Initial development release.
*   Core `NTSProject` framework.
*   Basic molecular networking and dual-network annotation.

