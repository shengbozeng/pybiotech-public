# Query Field Model

该模块（`pybiotech.classes.uniprot.https.uniprot.org.uniprot_query_field`）用 `pydantic.BaseModel` 描述 UniProt 查询字段的元数据，支持递归的 `QueryField` 结构与自动类型检查，方便上游根据 JSON/YAML 配置构造查询语句。

## 模块概览

| 名称 | 类型 | 描述 |
| --- | --- | --- |
| `QueryField` | `BaseModel` | 表示一个 UniProt 查询字段（term/label/type/regex 等），支持嵌套 `items`/`siblings` 用于构造分层 UI。 |
| `QueryField.model_rebuild` | class method | 重新构建前向引用，使嵌套的 `QueryField` 能互相引用而不会报 `ForwardRef` 错误。 |

---

### `QueryField`

**功能说明**

- 作为数据载体封装 UniProt 查询字段及其验证规则，字段定义与官方 JSON schema 保持对应。常见用法为 `QueryField(**entry)`，`entry` 来自 `xsdata`/JSON 解析结果。

**参数**

- `id` (`str`)：唯一 ID。
- `itemType` (`Literal["single","group","sibling_group"]`)：字段类别，决定如何拼接查询表达式。
- `label`, `term`, `dataType` 等字符串可选属性；`regex` 允许提供值校验模式；`siblings`/`items` 支持 `QueryField` 列表，实现递归树结构。

**返回类型**

- 返回 `QueryField` 实例；`pydantic` 会自动把嵌套结构转换为 `QueryField` 对象。

**异常**

- `pydantic.ValidationError`：数据类型/枚举不合法时会抛出；例如 `itemType` 非预设值或 `siblings` 不是可迭代。

**注意事项**

- ⚠️ `items` 与 `siblings` 可互相引用，务必在模块导入后调用 `model_rebuild()`（该模块在底部已执行）。
- ✅ 所有字段默认 `None`，便于原始 JSON 中缺失某些键时仍能成功创建。

**示例**

```python
from pybiotech.classes.uniprot.https.uniprot.org.uniprot_query_field import QueryField

field = QueryField(
    id="gene_name",
    itemType="single",
    term="protein_name",
    regex="^[A-Za-z0-9 ]+$",
    label="Protein name",
)
print(field.term)
```

---

### `QueryField.model_rebuild()`

**功能说明**

- 修复 `QueryField` 的前向引用，使 `siblings`/`items` 中的 `QueryField` 类型在运行时得以递归解析（这是 `pydantic` V2 的要求）。

**参数**

- 无；在模块导入阶段自动调用一次。

**返回类型**

- `None`；它直接修改类的前向引用元信息。

**异常**

- 不抛；但如果调用顺序错乱（例如先序列化再 `model_rebuild`）可能导致模型字段未正常绑定。

**注意事项**

- ℹ️ 已在模块末尾调用，使用者无需重复执行，除非在交互式环境中重新定义类。
