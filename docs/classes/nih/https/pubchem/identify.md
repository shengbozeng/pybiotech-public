# PubChem Identifier Helpers

`pybiotech.classes.nih.https.pubchem.identify` 提供面向 PubChem PUG REST 输入的标识符模型，负责清洗字符串、生成 URL 片段并按逗号处理 ID 列表。

## 模块概览

| 名称 | 类型 | 描述 |
| --- | --- | --- |
| `Identifier` | `BaseModel` | 表示单个标识符，支持 `/`/`&` 替换，附带 `type_hint` 用于后续逻辑。 |
| `Identifier.is_int` | method | 判断 `raw` 是否全数字。 |
| `Identifier.as_url_part` | method | 返回带 `quote` 转义的 URL 片段。 |
| `Identifiers` | `BaseModel` | 包装 `Identifier` 列表，支持 `from_string`（逗号分隔）与 `as_url_part`（生成拼接结果）。 |

---

### `Identifier`

**功能说明**

- 封装 PubChem 的单个 identifier，包含自动清洗（`/` 替 `.`，`&` 替 `%26`）与辅助方法用于 URL 拼接。

**参数**

- `raw` (`str`)：原始 ID；`model_validator` 当值含 `/` 或 `&` 时会自动转义。
- `type_hint` (`Optional[str]`)：可选，表示 `cid`/`smiles` 等（提升后续处理可读性）。

**返回类型**

- `Identifier` 实例。

**异常**

- `ValidationError`：若 `raw` 为空或非字符串。

**注意事项**

- ⚠️ `quote` 默认对 `/` 进行转义（与 `value.replace('/', '.')` 不冲突），确保生成合法 URL。

**示例**

```python
from pybiotech.classes.nih.https.pubchem.identify import Identifier

iden = Identifier(raw="DTP/NCI&something")
print(iden.as_url_part())  # DTP.NCI%26something
print(iden.is_int())  # False
```

---

### `Identifiers`

**功能说明**

- 表示逗号分隔的标识符序列，常用于 PubChem URL 的 `identifiers` 部分。

**参数**

- `items`：`Identifier` 列表；建议通过 `Identifiers.from_string` 来创建。

**返回类型**

- `Identifiers` 实例，`as_url_part` 会返回 `,` 分隔的个人 `Identifier` URL 片段。

**异常**

- 由底层 `Identifier` 的验证规则决定（例如非法字符会抛 `ValidationError`）。

**注意事项**

- ✅ `from_string` 会自动忽略空片段，确保 `items` 中只有合法字符串。

**示例**

```python
from pybiotech.classes.nih.https.pubchem.identify import Identifiers

ids = Identifiers.from_string("2244, 3672, DTP/NCI")
print(ids.as_url_part())  # 2244,3672,DTP.NCI
```
