# PubChem Input Spec

`pybiotech.classes.nih.https.pubchem.input_spec.InputSpec` 封装了 PUG REST 的 `<domain>/<namespace>/<identifiers>` 逻辑，提供自动猜测 `namespace`、支持 `parts` 自定义路径，并能反向解析 URL 片段。

## 模块概览

| 名称 | 类型 | 描述 |
| --- | --- | --- |
| `InputSpec` | `BaseModel` | 三段式（domain/namespace/identifiers）输入规格，包含自动规范化/校验与特殊 source 处理。 |
| `InputSpec.to_url_path` | method | 将字段拼接为合法 REST path，优先使用 `parts`（如 `["conformers", "..."]`）。 |
| `InputSpec.from_url_path` | class method | 反向解析 `domain/namespace/identifiers` 或任意自定义 path（返回 `parts`）。 |

---

### `InputSpec`

**功能说明**

- 描述 PubChem API 的输入片段，支持 `domain`（如 `compound`）、`namespace`（如 `cid`）、`identifiers`（字符串或列表），并针对 `sourceid` 等跳过 `/`、利用正则猜测命名空间。

**参数**

- `domain`：必须为 `DOMAIN_MAP` 中的值，`model_validator` 会将大小写统一为小写。
- `namespace`：若省略会从 `IDENTIFIER_GUESS` 中根据 `identifiers` 第一个值推断；`field_validator` 会校验命名空间是否在 `NAMESPACE_MAP` 中。
- `identifiers`：字符串或列表，支持逗号分隔；自动触发 `normalize_identifiers` 转为列表，并在 `namespace` 含 `source` 时把 `.` 转回 `/`。
- `parts`：优先级最高；若提供可处理像 `conformers/0000` 这类任意路径。

**返回类型**

- `InputSpec` 实例。

**异常**

- `ValueError`：域名非法、命名空间与域名不匹配、`identifiers` 类型错误。

**注意事项**

- ✅ `parts` 常用于 `conformers`、`annotations` 等路径；此时 `domain/namespace` 必须为空。
- ⚠️ `IDENTIFIER_GUESS` 只在通用输入时启用，不适用于 `sourceid` 或 `xref` 等特殊命名空间。

**示例**

```python
from pybiotech.classes.nih.https.pubchem.input_spec import InputSpec

spec = InputSpec(domain="compound", identifiers="C9H8O4")  # 自动推断为 formula
print(spec.to_url_path())  # compound/formula/C9H8O4

special = InputSpec(parts=["conformers", "000008C400000001"])
print(special.to_url_path())  # conformers/000008C400000001
```

---

### `InputSpec.to_url_path`

**功能说明**

- 将当前状态格式化为 URL 片段，优先使用 `parts`，其次拼接 `domain/namespace/identifiers`（会转义 `/` 与逗号）。

**返回类型**

- `str`：例如 `compound/cid/2244,962`。

**注意事项**

- ✅ 处理 `sourceid` 时自动把 `.` 替回 `/`，保证最终 URL 与 PubChem 要求一致。

**示例**

```python
spec = InputSpec(domain="substance", namespace="sourceid", identifiers="DTP/NCI")
print(spec.to_url_path())  # substance/sourceid/DTP.NCI
```

---

### `InputSpec.from_url_path(path: str)`

**功能说明**

- 读取 API path，若匹配标准三段式则返回相应字段；否则把整个分段写入 `parts`，便于解析例如 `annotations/foo`。

**参数**

- `path`：REST path（不含 base URL）。

**返回类型**

- `InputSpec`：`domain`/`namespace`/`identifiers` 或 `parts` 会被填充。

**异常**

- `ValueError`：当 path 无法解析为有效三段式时会 fallback 为 `parts`，不会直接抛。

**注意事项**

- ℹ️ `from_url_path` 会自动 decode `%XX` 字符并恢复 `.` 与 `/`。

**示例**

```python
spec = InputSpec.from_url_path("compound/cid/2244,962")
print(spec)
```
