# PubChem Namespace Compound

`pybiotech.classes.nih.https.pubchem.section.nscompound.NS` 是 `NamespaceModel` 的具体实现，用于链接 compound/substance 等 domain 与其 namespace，同时支持前缀匹配（如 `fastidentity/xxxx`）和值替代。

## 模块概览

| 名称 | 类型 | 描述 |
| --- | --- | --- |
| `NS` | class | 继承自 `NamespaceModel`，实现 `validate_namespace_and_domain` 并允许 namespace 以 `allowed/` 前缀开头。 |

---

### `NS`

**功能说明**

- 校验 `value` 是否在 `DOMAIN_NAMESPACE_MAP[domain]` 范围内，或以其中某个 namespace 作为前缀（例如 `fastidentity/xxx`）后再生成可用的 `value`。

**参数**

- 继承 `NamespaceModel` 的 `value`（namespace 字符串）与 `domain`。

**返回类型**

- `NS` 实例；`as_url_part()` 直接返回 `value`。

**异常**

- `ValueError`：当 namespace 未被允许且不以任意 `allowed/` 前缀开头时抛出。

**注意事项**

- ✅ 支持 `fastidentity`、`smiles` 等组合形式；只要 `value` 以 `allowed` 之一 + `/` 为前缀即可通过验证。

**示例**

```python
from pybiotech.classes.nih.https.pubchem.section.domain import ENUMDomain
from pybiotech.classes.nih.https.pubchem.section.nscompound import NS

ns = NS(value="cid", domain=ENUMDomain.compound)
print(ns.as_url_part())
```
