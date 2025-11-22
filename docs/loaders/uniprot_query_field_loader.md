# UniProt Query Field Loader

该模块封装了 UniProt 查询字段的读取、验证与构造，依赖 `QueryField` dataclass 表达每个字段定义，适合基于 JSON 配置快速构造查询字符串。

## 模块概览

| 名称 | 类型 | 描述 |
| --- | --- | --- |
| `UniProtQueryFieldLoader` | class | 静态方法集合，围绕 `QueryField` 的加载、格式验证与查找构建工具。 |
| `UniProtQueryFieldLoader.load_query_fields` | static method | 从文件/字符串/字节载入 QueryField 定义，返回解析后的列表。 |
| `UniProtQueryFieldLoader.build_query` | static method | 用 `QueryField` 和用户输入构建合法的 UniProt 查询片段（term:value）。 |
| `UniProtQueryFieldLoader.find_field_by_term` | static method | 在字段树（包括 `siblings`/`items`）中递归查找某个 `term`。 |

---

### `UniProtQueryFieldLoader.load_query_fields(json_input: Union[Path, str, bytes]) -> List[QueryField]`

**功能说明**

- 把 JSON 输入解析为 `QueryField` 实例列表，适合集群运行时把官方 JSON 配置转为 Python 对象。

**参数**

- `json_input`：可以是 `Path`（必须指向已存在的文件）、JSON 字符串或字节；文件会以 `utf-8` 打开读取。

**返回类型**

- `List[QueryField]`：顺序保持原 JSON；每个元素都通过 `QueryField(**entry)` 初始化。

**异常**

- `FileNotFoundError`：文件路径不存在。
- `TypeError`：`Path` 不是文件、或输入类型不在 `Path/str/bytes` 里。
- `json.JSONDecodeError`：字符串/文件内不是合法 JSON。

**注意事项**

- 文件读取时使用 `encoding="utf-8"`，JSON 结构必须是列表；出现额外键会传给 `QueryField`，避免漏掉默认值。

**示例**

```python
from pathlib import Path
from pybiotech.loaders.uniprot_query_field_loader import UniProtQueryFieldLoader

fields = UniProtQueryFieldLoader.load_query_fields(Path("config/uniprot_query_fields.json"))
```

---

### `UniProtQueryFieldLoader.build_query(field: QueryField, value: str) -> str`

**功能说明**

- 根据 `QueryField` 的定义验证 `value`（挡住正则）、然后输出标准的 `term:value` 或其 `siblings` 的 `term`，保证生成的查询字符串可以直接贴给 UniProt HTTP API。

**参数**

- `field`：已经加载好的 `QueryField`，`regex`/`itemType`/`term`/`siblings` 等字段决定行为。
- `value`：要插入的查询值，会被 `str()` 强制转为字符串再判断正则。

**返回类型**

- `str`：合法的 UniProt 查询片段，例如 `gene:TP53`。

**异常**

- `ValueError`：`value` 不满足 `field.regex`，会附带 `field.example`（如有）作为提示。
- `NotImplementedError`：碰到未知的 `itemType`（例如不是 `single` 且没有 `siblings`）。

**注意事项**

- `itemType == "single"` 时直接使用 `field.term`；`"sibling_group"` 会取第一个 `siblings[0]` 的 term。其他类型目前不支持。

**示例**

```python
field = fields[0]
query_segment = UniProtQueryFieldLoader.build_query(field, "P04637")
```

---

### `UniProtQueryFieldLoader.find_field_by_term(fields: List[QueryField], term: str) -> Optional[QueryField]`

**功能说明**

- 在字段树中依据 `term` 递归查找 `QueryField`，允许在 `siblings` 与 `items` 子树里查找，返回第一条匹配。

**参数**

- `fields`：从 `load_query_fields` 返回的列表或任意一组 `QueryField`。
- `term`：待查找的字段名称，会与 `field.term` 及其 `siblings`/`items` 的 `term` 比较。

**返回类型**

- `Optional[QueryField]`：找到了返回对象，否则返回 `None`。

**异常**

- 无，找不到只返回 `None`，不抛。

**注意事项**

- 使用递归，先检查当前层，再遍历 `siblings` & `items`；若结构较深，可能会多层调用，但 `QueryField` 树通常较浅。

**示例**

```python
field = UniProtQueryFieldLoader.find_field_by_term(fields, "protein_name")
if field:
    print(field.description)
```
