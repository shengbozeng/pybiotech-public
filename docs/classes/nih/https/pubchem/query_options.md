# PubChem Query Options

`pybiotech.classes.nih.https.pubchem.query_options.QueryOptions` 处理 PUG REST 的 query string，支持官方常见参数、类型校验、互斥规则、JSONP 约束，并提供正/反向解析。

## 模块概览

| 名称 | 类型 | 描述 |
| --- | --- | --- |
| `QueryOptions` | class | 管理 query 参数字典，内置 `_OPTIONS_SCHEMA`、`set_context`、`contextual_validate`。 |
| `QueryOptions.set_context` | method | 注入 `operation`/`output_format` 上下文以供后续校验。 |
| `QueryOptions.contextual_validate` | method | 执行参数互斥、依赖与 allowable-values 规则。 |
| `QueryOptions.to_query_string` | method | 生成符合 RFC 3986 的 query string。 |
| `QueryOptions.from_query_string` | class method | 从 URL query 反解为 `QueryOptions` 实例。 |

---

### `QueryOptions`

**功能说明**

- 接受任意 PubChem query 参数，尝试用 `_OPTIONS_SCHEMA` 做基本类型校验（字符串、布尔、整数、枚举列表），同时保留自定义字段，最后生成 `?key=value` 形式的 query string。

**参数**

- 任意 `key=value` 组合：构造时会通过 `__setitem__` 调用 `_validate`，如果 key 在 `_OPTIONS_SCHEMA` 中则会做类型/值校验（例如 `record_type` 只能为 `2d/3d`）。

**返回类型**

- `QueryOptions` 实例，内部 `params` 字典可通过 `items()` 访问。

**异常**

- `ValueError`：参数值不在允许范围（例如 `record_type="4d"`）或 `callback` 与 format 不匹配等。

**注意事项**

- ✅ `_OPTIONS_SCHEMA` 列出了 30+ 常见选项，包含布尔、枚举、整型；对于未声明字段仍然允许（便于扩展）。
- ⚠️ `__setitem__` 会调用 `_allow_special`，比如 `image_size=320x240` 会被接受。

**示例**

```python
opts = QueryOptions(record_type="3d", image_size="320x240", MaxRecords=50)
```

---

### `QueryOptions.set_context(operation: Optional[str], output_format: Optional[str])`

**功能说明**

- 记录当前 request 的 operation/output_format，以便 `contextual_validate` 中进行业务规则校验。

**参数**

- `operation`：如 `property`/`record`。
- `output_format`：如 `JSON`/`PNG`。

**注意事项**

- ✅ 推荐在 `PubChemRequest` 初始化时调用，确保 `contextual_validate` 有足够上下文。
---

### `QueryOptions.contextual_validate()`

**功能说明**

- 校验参数之间的互斥或依赖关系（如 `callback` 只能与 `JSONP` 一起出现、`PN` 仅适用于 `record`/`conformers`、`sid=listkey` 必须有 `listkey` 等）。

**返回类型**

- `None`；校验失败会抛 `ValueError`。

**注意事项**

- ⚠️ 一定在构建 URL 之前调用（`PubChemRequest.validate` 会自动调用）；若 catch `ValueError` 需向调用者提示具体参数冲突。

**示例**

```python
opts.set_context(operation="record", output_format="PNG")
opts.contextual_validate()
```

---

### `QueryOptions.to_query_string()`

**功能说明**

- 把参数字典编码为 `?key=value&...`，默认使用 `quote_plus`，确保空格与特殊字符正确编码。

**返回类型**

- `str`：包括前导 `?`，空参数字典返回空字符串。

**示例**

```python
print(QueryOptions(record_type="3d").to_query_string())  # ?record_type=3d
```

---

### `QueryOptions.from_query_string(qs: str)`

**功能说明**

- 反向解析 query string（可带 `?`），合并重复字段（例如 `sid=1&sid=2` 会转成 `sid=1,2`）。

**返回类型**

- 新的 `QueryOptions` 实例。

**示例**

```python
parsed = QueryOptions.from_query_string("?record_type=3d&image_size=320x240")
```
