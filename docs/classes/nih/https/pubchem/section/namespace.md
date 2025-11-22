# PubChem Section Namespace

`pybiotech.classes.nih.https.pubchem.section.namespace` 定义了 PubChem domain 关联的 namespace 枚举与验证基础类，用于确保 namespace（如 `cid`/`smiles`）在该 domain 内合法。

## 模块概览

| 名称 | 类型 | 描述 |
| --- | --- | --- |
| `ENUNameSpace` | `Enum` | PubChem 支持的 namespace 集合（`cid`、`smiles`、`inchikey`、`listkey` 等）。 |
| `DOMAIN_NAMESPACE_MAP` | `dict` | 映射 `ENUMDomain` 到可接受 `ENUNameSpace` 列表。 |
| `NamespaceModel` | `BaseModel` | 抽象基类，提供 `value` + `domain` 字段，并定义 `validate_namespace_and_domain` 与 `as_url_part()`。 |

---

### `ENUNameSpace`

**功能说明**

- 枚举 PubChem `namespace`，便于在请求构造时进行简单枚举校验。

**示例**

```python
from pybiotech.classes.nih.https.pubchem.section.namespace import ENUNameSpace

print(ENUNameSpace.cid.value)
```

---

### `DOMAIN_NAMESPACE_MAP`

**功能说明**

- 将 `ENUMDomain`（例如 `compound`）映射到合法的 `ENUNameSpace` 列表；用于 `NS` 子类的 `validate_namespace_and_domain`。

**示例**

```python
from pybiotech.classes.nih.https.pubchem.section.namespace import DOMAIN_NAMESPACE_MAP

print(DOMAIN_NAMESPACE_MAP.keys())
```

---

### `NamespaceModel`

**功能说明**

- 抽象基类，要求子类实现 `validate_namespace_and_domain`（`@abstractmethod`），并提供 `as_url_part()` 方便拼接 URL 片段。

**参数**

- `value`：`str` namespace 字符串。
- `domain`：`ENUMDomain`。

**返回类型**

- 继承自 `NamespaceModel` 的具体实现（例如 `NS`）。

**注意事项**

- ℹ️ `validate_namespace_and_domain` 的实现通常会调用 `DOMAIN_NAMESPACE_MAP` 验证 `value` 是否允许。

**示例**

```python
from pybiotech.classes.nih.https.pubchem.section.namespace import ENUNameSpace, DOMAIN_NAMESPACE_MAP
from pybiotech.classes.nih.https.pubchem.section.domain import ENUMDomain

print(DOMAIN_NAMESPACE_MAP[ENUMDomain.compound])
```
