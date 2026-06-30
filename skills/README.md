# Hermes Skills Library

小賀的「程序知識」庫。Skills 是 Workflows 的補充，封裝的是**決策邏輯**和**診斷推理**，
而非執行步驟（執行步驟在 `workflows/` 的 Snakefile 裡）。

## 架構說明

```
skills/                          ← 程序知識（此目錄）
├── diagnose-snakemake-failure.md  ← 如何診斷錯誤
├── choose-variant-caller.md       ← 如何選工具/配置
├── interpret-qc-metrics.md        ← 如何解讀結果
└── ...

workflows/                       ← 執行步驟（Snakemake DAG）
├── wgs_snp/Snakefile
├── rnaseq/Snakefile
└── ...
```

## 現有 Skills

| Skill | 觸發時機 |
|-------|---------|
| [intake-qc-checklist](./intake-qc-checklist.md) | **用戶描述新分析需求時（最先觸發）**：資料完整性、樣本設計、各 workflow 特定紅線 |
| [diagnose-snakemake-failure](./diagnose-snakemake-failure.md) | Workflow 失敗、crash、`PermissionError`、`MissingOutputException` |
| [choose-variant-caller](./choose-variant-caller.md) | 用戶描述 SNP/INDEL 分析需求 |
| [interpret-qc-metrics](./interpret-qc-metrics.md) | 用戶問「QC 結果如何？」「可以繼續嗎？」 |

## 如何使用

小賀在對話中遇到對應情境時，主動參考對應的 skill 文件。
Skills 不是腳本，是「大腦的決策流程」——告訴小賀**怎麼想**，而非**做什麼**。

## 如何新增 Skill

1. 在此目錄新增 `<skill-name>.md`
2. 用條件句描述觸發時機（「當用戶問 X 時」）
3. 包含決策樹、參考數字、常見錯誤的處理方式
4. 更新此 README 的表格

## 與 ClawBio 的關係

ClawBio（OpenClaw 框架）的 skills 格式為 SKILL.md + frontmatter，
本目錄採用更輕量的純 Markdown 格式，直接適配 Hermes/Claude 的對話環境，
無需 OpenClaw agent runtime。

若未來引入 ClawBio skills，可將其邏輯（決策樹、最佳實踐）改寫為此格式，
保留知識但不依賴其執行引擎。
