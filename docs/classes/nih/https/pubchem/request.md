# PubChem Request Builder

`pybiotech.classes.nih.https.pubchem.request.PubChemRequest` 封装了构建 PubChem PUG REST URL、Headers、Body 及响应解释的全流程，聚合 `InputSpec`、`OperationSpec`、`OutputSpec` 与 `QueryOptions`。

## 模块概览

| 名称 | 类型 | 描述 |
| --- | --- | --- |
| `PubChemRequest` | class | 核心构造器，接受输入/操作/输出/QueryOptions，自动注入上下文并提供 `build_url`、`build_headers`、`build_body`、`process_response` 等辅助方法。 |
| `validate` | method | 校验 `QueryOptions` 上下文、`operation`/`output` 的组合合法性、JSONP callback、`property` 操作必须带 tags 等。 |
| `build_url` | method | 先调用 `validate()`，再拼接 `input_spec`、`operation_spec`、`output_spec` 与 query string。 |
| `build_headers` | method | 根据输出格式映射 `Accept`，POST 请求附加 `Content-Type`。 |
| `build_body` | method | 必要时返回 URL 编码的 POST body；否则 `None`。 |
| `process_response` | method | 简单映射 HTTP 状态码到描述，并返回结构化字典。 |
| `from_url` | class method | 解析完整 REST URL，逐步重建 `InputSpec`/`OperationSpec`/`OutputSpec`/`QueryOptions`。 |

---

### `PubChemRequest`

**功能说明**

- 把输入/操作/输出/选项组合成合法 request，内部会调用 `query_options.set_context(operation, output_format)` 并在构造时缓存 `use_post`/`post_body`。

**参数**

- `input_spec`：`InputSpec`，负责 `domain/namespace/identifiers`。
- `operation_spec`：`OperationSpec`，定义 `<operation>`；可为 `None`（默认 `record`）。
- `output_spec`：`OutputSpec`，控制输出格式/JSONP callback。
- `query_options`：可选，默认新建 `QueryOptions()`。
- `use_post`、`post_body`：控制是否以 POST 发送。

**返回类型**

- `PubChemRequest` 实例。

**注意事项**

- ✅ 构造时即调用 `query_options.set_context` 使 `contextual_validate` 能判断 `operation` 与 `output_format`。
---

### `validate()`

**功能说明**

- 校验所有上下文组合，包括：
  - `query_options.contextual_validate()`；
  - `property` 必须带 `tags`；
  - `PNG/ TXT` 与特定 `operation` 互斥；
  - `JSONP` 必须提供 callback。

**异常**

- `ValueError`：任何违规组合都会抛出，提示具体规则（例如 `PNG` 只能搭配 `record/conformers`）。

**注意事项**

- ⚠️ `build_url` 里会隐式调用 `validate()`，无需外部重复调用。

---

### `build_url(base_url: str = "https://pubchem.ncbi.nlm.nih.gov/rest/pug")`

**功能说明**

- 在合法参数下拼接 `base_url + input + operation + output + querystring`（`query_options.to_query_string()`），支持 `operation_spec` 为 `None` 的情况（仅 `input`）。

**返回类型**

- 完整 REST URL（`str`）。

**示例**

```python
req = PubChemRequest(input_spec, operation_spec, output_spec)
print(req.build_url())
```

---

### `build_headers()`

**功能说明**

- 根据 `output_spec.output_format` 映射官方 `Accept` 头（JSON → `application/json`），若 `use_post` 设置为 `True`，还会加入 `Content-Type: application/x-www-form-urlencoded`。

**返回类型**

- `dict`，可直接传给 `requests`。

---

### `build_body()`

**功能说明**

- 若 `use_post` 为 `True` 且提供 `post_body` 字典，会用 `urlencode`（`quote_plus`）返回字符串；否则返回 `None`。

**注意事项**

- ✅ 确保 `post_body` 的键/值已经安全格式化。

---

### `process_response(status_code: int, response_text: str)`

**功能说明**

- 简单映射常见 HTTP 状态码到友好描述，返回 `dict`（包含 `raw` 文本），便于上层做统一日志。

**返回类型**

- `dict`：`{"status_code": ..., "status": "...", "raw": ...}`。

---

### `from_url(url: str)`

**功能说明**

- 解析完整 REST URL：通过正则剥离 base_url，再依次使用 `InputSpec.from_url_path`、`OperationSpec.from_url_path`、`OutputSpec.from_url_path`、`QueryOptions.from_query_string` 重建对象。

**异常**

- `ValueError`：当 URL 不符合 `https://.../rest/pug/...` 或无法识别其中某个片段时抛出。

**示例**

```python
req = PubChemRequest.from_url("https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/2244/JSON")
```
