# PubChem Output Spec

`pybiotech.classes.nih.https.pubchem.output_spec.OutputSpec` 负责表示 `<output specification>`（如 `JSON`、`SDF`、`JSONP?callback=...`），包含合法格式校验、JSONP callback 处理，以及与 URL 片段的互转。

## 模块概览

| 名称 | 类型 | 描述 |
| --- | --- | --- |
| `OutputSpec` | `BaseModel` | 保存 `output_format` 和 `callback`，并在 `model_validator` 中统一大写、校验并清理不合法组合。 |
| `OutputSpec.OUTPUT_FORMATS` | `ClassVar[list]` | 支持的格式列表（`"XML"`、`"JSON"`、`"SDF"`、`"PNG"` 等）。 |
| `OutputSpec.to_url_path` | method | 输出合法的 URL 片段，如 `JSONP?callback=my_cb`。 |
| `OutputSpec.from_url_path` | class method | 从 URL 片段解析出格式 + callback。 |

---

### `OutputSpec`

**功能说明**

- 统一表示 PubChem 输出段，确保 `output_format` 合法（在 `OUTPUT_FORMATS` 内）并强制把 `callback` 仅在 `JSONP` 时保留。

**参数**

- `output_format`：字符串（大小写不敏感）；`model_validator` 会转为大写并校验列表。
- `callback`：JSONP 回调名称；非 `JSONP` 格式时会被自动置 `None`。

**返回类型**

- `OutputSpec` 实例。

**异常**

- `ValueError`：格式不在白名单内（例如 `output_format="ZIP"` 会报错）。

**注意事项**

- ✅ `model_validator` 内部默认只接受官方支持的枚举格式，并自动清除 `callback`。

**示例**

```python
from pybiotech.classes.nih.https.pubchem.output_spec import OutputSpec

spec = OutputSpec(output_format="json")
print(spec.output_format)  # JSON
```

---

### `OutputSpec.to_url_path()`

**功能说明**

- 把 `OutputSpec` 转为 REST path；若是 `JSONP` 且提供 callback，会附加 `?callback=`。其他格式直接返回大写名称。

**返回类型**

- `str`：例如 `JSON`, `JSONP?callback=my_cb`, `SDF`，空字符串表示未指定格式。

**注意事项**

- 🧭 `quote` 只在 `callback` 上使用，防止特殊字符破坏 URL。

**示例**

```python
spec = OutputSpec(output_format="jsonp", callback="cb")
print(spec.to_url_path())  # JSONP?callback=cb
```

---

### `OutputSpec.from_url_path(path: str)`

**功能说明**

- 反向解析 `JSONP?callback=foo` 或简单 `CSV` 路径，带参数的格式会把 `callback` 提取出来。

**参数**

- `path`：`output` 段（可含 `?`）。

**返回类型**

- `OutputSpec`；如果 `path` 为空则返回默认实例。

**注意事项**

- ⚠️ 若路径在 `?` 之后还带其他查询参数，当前实现仅识别 `callback`。

**示例**

```python
spec = OutputSpec.from_url_path("JSONP?callback=cb")
```
