# PubChem Section Domain

`pybiotech.classes.nih.https.pubchem.section.domain` 明确了 PubChem `domain` 与其 URL 表示，提供 `ENUMDomain` 枚举和 `DomainModel` 供上游统一接口使用。

## 模块概览

| 名称 | 类型 | 描述 |
| --- | --- | --- |
| `ENUMDomain` | `Enum` | 列出 PubChem 支持的 domain（`compound`、`substance`、`assay`、`gene` 等）。 |
| `DomainModel` | `BaseModel` | 包装 `ENUMDomain`，提供 `as_url_part()` 生成 path 片段并通过 `field_validator` 确保合法值。 |

---

### `ENUMDomain`

**功能说明**

- 将常用 domain 封装为枚举，用于请求构造与校验。

**返回类型**

- 枚举值，如 `ENUMDomain.compound`。

**注意事项**

- ✅ `DomainModel` 会把输入（如 `"Compound"`）规范化为 `compound`。

**示例**

```python
from pybiotech.classes.nih.https.pubchem.section.domain import ENUMDomain

print(ENUMDomain.compound.value)  # compound
```

---

### `DomainModel`

**功能说明**

- 以 `ENUMDomain` 为类型包装 `value`，`field_validator` 会把传入的字符串转成枚举；`as_url_part()` 输出 `value.value` 供 URL 直拼。

**参数**

- `value`：`str` 或 `ENUMDomain`，会统一通过 validator 验证。

**返回类型**

- `DomainModel` 实例。

**异常**

- `ValueError`：传入非法 domain（如 `foo`）时抛出。

**示例**

```python
from pybiotech.classes.nih.https.pubchem.section.domain import DomainModel

model = DomainModel(value="Compound")
print(model.as_url_part())  # compound
```
