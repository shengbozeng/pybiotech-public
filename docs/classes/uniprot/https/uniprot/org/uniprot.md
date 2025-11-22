# UniProt Schema

该模块 (`pybiotech.classes.uniprot.https.uniprot.org.uniprot`) 是 `xsdata` 从 UniProt XML 架构生成的 dataclass/enum 集合，覆盖 Entry、Protein、Sequence 等主要元素，所有属性都已标注类型与 namespace，用于高保真地内存表示 UniProt XML 解析结果。

## 模块概览

| 名称 | 类型 | 描述 |
| --- | --- | --- |
| `__NAMESPACE__` | `str` | 模块级常量，统一为 `http://uniprot.org/uniprot`；所有 dataclass 都引用该 namespace。 |
| `EntryDataset` | `Enum` | 表示 `entry/@dataset` 的来源——`Swiss-Prot`/`TrEMBL`。 |
| `Entry` | `@dataclass` | UniProtKB 的根结构，包含 accession、name、protein、organism、comment、sequence 等字段。 |
| `ProteinType` | `@dataclass` | 包含 `recommendedName`、`alternativeName`、`submittedName` 等子结构，等价于扁平化的 DE 行。 |
| `ProteinExistenceType` | `@dataclass` | 对应 PE 行，结合 `ProteinExistenceTypeType` 枚举指出证据等级。 |
| `SequenceType1` | `@dataclass` | 记录序列字符串 + `length`/`mass`/`checksum`/`modified` 等属性。 |
| `Controlled enums` | `Enum` | 一组由 XSD 生成的枚举（`CitationTypeType`、`CommentTypeType`、`FeatureTypeType`、`PositionTypeStatus`、`StatusTypeStatus` 等），用于保证字段值在预定义词汇之内。 |

---

### `EntryDataset`

**功能说明**

- 提供 `entry` 标签的 `dataset` 属性，对应两个主库：`Swiss-Prot` 与 `TrEMBL`，可用于区分手工审核与自动注释记录。

**返回类型**

- `EntryDataset` 枚举：`.SWISS_PROT` 或 `.TR_EMBL`。

**注意事项**

- ℹ️ `xsdata` 会将 XML 属性转换为 `EntryDataset`，无需手动映射字符串；直接访问 `entry.dataset` 即可。

**示例**

```python
from pybiotech.classes.uniprot.https.uniprot.org.uniprot import EntryDataset

print(EntryDataset.SWISS_PROT.value)  # Swiss-Prot
```

---

### `Entry`

**功能说明**

- 是 UniProtKB 的核心结构，表示一个 `entry` 标签，封装 `accession`/`name` 列表、`ProteinType`、`organism`、`reference`、`comment`、`sequence` 等信息，适合用于数据库导入或高保真序列解析。

**参数**

- `accession`（`list[str]`）：至少一个 accession。
- `name`（`list[str]`）：可包含主名与别名。
- `protein`（`ProteinType`）：详述 protein 名称系。
- `organism`/`organism_host`/`gene`：均为各自 dataclass 列表，详见模块。
- `comment`、`feature`、`db_reference`、`evidence`：可选列表，`CommentType`、`FeatureType` 等字段均采用 `Enum`/`EvidencedStringType` 增强。
- `sequence`（`SequenceType1`）：必需的序列数据（含 checksum、length 等）。
- `dataset`（`EntryDataset`）、`created`/`modified`（`XmlDate`）、`version`：来自 XML 属性。

**返回类型**

- `Entry` 实例。

**异常**

- `TypeError`/`ValueError`：若缺少 `accession` 或 `sequence`，`xsdata` 在构造时会抛出。

**注意事项**

- ⚠️ `entry` 可能含多个 `comment`/`feature`，使用 `entry.comment` 迭代可获取每个 `CommentType`（该类会附带 `type`/`text`）。
- ✅ 构造实例常见方式是用 `xsdata` 直接从 XML 反序列化：`XmlParser(...).from_string(xml_str, Entry)`。

**示例**

```python
from xsdata.formats.dataclass.parsers import XmlParser
from xsdata.formats.dataclass.context import XmlContext
from pybiotech.classes.uniprot.https.uniprot.org.uniprot import Entry

parser = XmlParser(context=XmlContext())
entry = parser.from_string(xml_snippet, Entry)
print(entry.accession[0], entry.sequence.length)
```

---

### `ProteinType`

**功能说明**

- 表达 DE 行的多种名称（`recommendedName`、`alternativeName`、`submittedName`、`allergenName` 等），每个子成员都支持 `EvidencedStringType` 以保留证据来源。

**参数**

- `recommended_name`：包含 `fullName`、`shortName`、`ecNumber` 等字段的 `ProteinType.RecommendedName`。
- `alternative_name`/`submitted_name`：可选列表，每项记录 `fullName` + `ecNumber`/`synonym`。
- `allergen_name`、`biotech_name` 等为 `EvidencedStringType`，附带 `evidence` ID 列表。

**返回类型**

- `ProteinType` 实例。

**注意事项**

- ℹ️ `EvidencedStringType` 的 `evidence` 属性可用于追踪注释来源。
- ✅ `ProteinType` 中的 `recommendedName.full_name.value` 即可取得主名称。

**示例**

```python
protein = ProteinType(
    recommended_name=ProteinType.RecommendedName(
        full_name=EvidencedStringType(value="Hemoglobin subunit beta")
    )
)
```

---

### `SequenceType1`

**功能说明**

- 保存完整序列字符串与物理属性（`length`、`mass`、`checksum`、`modified`、`version`）以及 `fragment`/`precursor` 标记，等价于 flat-file 的 SEQ 行。

**参数**

- `value`：完整的氨基酸序列（`str`）。
- `length`/`mass`/`checksum`：标量属性；都标记为必需。
- `modified`：`XmlDate`，表示该版本提交时间。
- `fragment`：若序列为片段，`SequenceTypeFragment` 会携带 `status`。

**返回类型**

- `SequenceType1` 实例。

**注意事项**

- ⚠️ `value` 可能包含换行；使用 `''.join(entry.sequence.value.split())` 提取纯序列。
- ✅ `checksum` 可用于与 PubMed/UniProt REST 验证一致性。

**示例**

```python
seq = SequenceType1(
    value="MVHLTPEEKSAVTALWGKVNVDEVGGEALGRLLVVYPWTQR",
    length=40,
    mass=10455,
    checksum="F8F4F0DB",
    modified=XmlDate(2025, 1, 1),
    version=1,
)
```
