# PubChem Operation Spec

`pybiotech.classes.nih.https.pubchem.opera_spec.OperationSpec` 用于描述 REST URL 中的 `<operation specification>`（例如 `property/MolecularFormula`、`targets/GeneID`），具备 domain 校验、tag 自动拆分以及对 multi-level operation 的序列化/反序列化能力。

## 模块概览

| 名称 | 类型 | 描述 |
| --- | --- | --- |
| `OperationSpec` | `BaseModel` | 记录 `operation`、可选 `tags`、可选 `domain` 并在 `model_validator` 中保证合法组合。 |
| `OperationSpec.OPERATION_MAP` | `ClassVar[Dict[str, List[str]]]` | 每个 domain 允许的 operation（如 `compound` 允许 `record`、`property`、`xrefs` 等）。 |
| `OperationSpec.MULTI_LEVEL_OP` | `ClassVar[List[str]]` | `property`、`xrefs`、`targets` 等多级 operation，格式：`operation/tag1,tag2`。 |
| `OperationSpec.to_url_path` | method | 根据 operation + tags 生成 path，自动对 tag 列表 `quote`。 |
| `OperationSpec.from_url_path` | class method | 将格式如 `property/MolecularFormula,InChIKey` 解析为 `OperationSpec`。 |
| `OperationSpec.validate` | method | 显式触发 `pydantic` 校验（通常不必调用）。 |

---

### `OperationSpec`

**功能说明**

- 表示 `<operation>` 段，包括支持单层 `record`、`synonyms`，以及多层 `property/MolecularFormula,InChIKey`，并使用 `domain` 进行上下文校验。

**参数**

- `operation`：必填（如 `property`、`targets/GeneID`），会被转义以构成 URL。
- `tags`：可为字符串（`"MolecularFormula,InChIKey"`）或列表；`model_validator` 会拆分并清洗为空字符串的项。
- `domain`：辅助校验的上下文，`normalize_domain` 会转小写再校验 `OPERATION_MAP`。

**返回类型**

- `OperationSpec` 实例。

**异常**

- `ValueError`：当 `operation` 不在 `OPERATION_MAP[domain]` 中，如 `operation="png"` 但 domain=`compound` 时。
- `ValueError`：`operation="property"` 但 `tags` 为空。

**注意事项**

- ✅ `tags` 支持列表/字符串，`auto_split_tags` 会在构造后把逗号切分。
- ⚠️ `domain` 为 `assay` 时可构造 `targets/GeneID`，`OperationSpec` 也允许 `operation` 自带 `/`（例如 `targets/ProteinGI`）。

**示例**

```python
from pybiotech.classes.nih.https.pubchem.opera_spec import OperationSpec

spec = OperationSpec(operation="property", tags=["MolecularFormula", "InChIKey"], domain="compound")
```

---

### `OperationSpec.to_url_path()`

**功能说明**

- 生成拼接到 `input` 后的 path 部分，`property`/`xrefs`/`targets` 会把 tag 列表 `join` 并 `quote`。

**返回类型**

- `str`：例如 `property/MolecularFormula,InChIKey` 或 `record`。

**注意事项**

- 🧭 若 `operation` 本身已经包含 `/`（如 `targets/GeneID`），则直接返回该值。

**示例**

```python
print(OperationSpec(operation="record", domain="compound").to_url_path())
```

---

### `OperationSpec.from_url_path(path: str, domain: Optional[str] = None)`

**功能说明**

- 反向恢复 `OperationSpec`，支持 `tag` 字符串的拆分（逗号/空白）并根据 `domain` 校验合法性。

**参数**

- `path`：如 `property/MolecularFormula,InChIKey`。
- `domain`：可选；用于 `OPERATION_MAP` 校验。

**返回类型**

- `OperationSpec`。

**异常**

- `ValueError`：`path` 为空、`domain` 不支持指定操作等会抛出。

**示例**

```python
spec = OperationSpec.from_url_path("property/MolecularFormula,InChIKey", domain="compound")
```
