# 如何给分析代码获取一个 DOI（Zenodo 免费流程）

当前版本：**v1.1.0**（2026-08-18；在 v1.0.0 / 10.5281/zenodo.21947642 基础上新增
GSE280220 独立验证、ODE 汇总、Enrichr 全量表、固定种子的 MR-PRESSO R5 输出）。
仓库根目录已包含 `.zenodo.json`（标题、作者、ORCID、关键词、许可证、版本、关联
v1.0.0 DOI 等元数据），GitHub 集成会自动读取。

## 方法一：GitHub + Zenodo 自动关联（推荐）

1. 把本文件夹（`analysis_repository`）上传为一个 GitHub 仓库：`https://github.com/XuZhengdong2026/pyoderma-gangrenosum-mr`。
2. 打开 https://zenodo.org ，用 GitHub 账号登录（首次需在 Zenodo 页面点击 "GitHub" 并授权）。
3. 在 Zenodo 的 GitHub 页面找到该仓库，把开关拨到 "On"（开启自动归档）。
4. 把 `.zenodo.json` 一并提交到仓库根目录（本文件夹已包含）。
5. 回到 GitHub 仓库，创建一次 Release（例如 `v1.1.0`，tag 与版本号一致）。
   Zenodo 会自动为新版本生成一个 DOI；如果该仓库之前已关联 v1.0.0 的存档，
   新 Release 会作为 v1.0.0 的新版本（new version）归档。
6. 把新 DOI 填进稿件 Data Availability 声明即可。

注意：每次发新 Release 会生成新版本 DOI；稿件中写明版本号（v1.1.0）和所用 DOI。

## 方法二：直接上传 Zenodo（不需要 GitHub）

1. 打开 https://zenodo.org ，注册/登录。
2. 点击 "New upload"，上传本文件夹打包的 zip
   （`analysis_repository_v1.1.0_20260818.zip`，位于
   `F:\坏疽性脓皮病\outputs\`）。
3. 按 `.zenodo.json` 内容填写标题、作者（含 ORCID）、描述、关键词，
   资源类型选 **Software**，版本填 **v1.1.0**，许可证选 **MIT**，
   Access 选 **Open Access**。
4. 在 Related/alternate identifiers 中添加
   `10.5281/zenodo.21947642`，关系选 **isPreviousVersionOf**。
5. 保存后点击 "Publish"，页面会给出新 DOI。

## 稿件中建议的写法

> Analysis code for instrument selection, MR, mediation, specificity, transcriptomic (including the independent GSE280220 validation), network and reproducibility analyses is available at https://github.com/XuZhengdong2026/pyoderma-gangrenosum-mr (version v1.1.0; DOI: 10.5281/zenodo.21991327).

把版本号和 DOI 替换成实际值后，填入两份稿件的 Data Availability 声明即可。
