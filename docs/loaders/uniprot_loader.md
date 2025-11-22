# UniProt Loader

该模块围绕 `pybiotech.classes.uniprot.https.uniprot.org.uniprot.Entry` 提供 XML 输入的流式迭代解析能力，支持文件、字符串或字节序列，并能自动修正 namespace。

## 模块概览

| 名称 | 类型 | 描述 |
| --- | --- | --- |
| `UniProtLoader` | class | 通过 `xsdata` 的 `XmlParser` 将 UniProt XML 转换为 `Entry` dataclass；构造过程中会自动识别定义 namespace。 |
| `UniProtLoader.iterate_entries` | instance method | 按顺序从输入中逐条生成 `Entry` 对象；文件版使用 `lxml.etree.iterparse`，字符串版重写 namespace 后逐一解析。 |
| `UniProtLoader.detect_namespace` | static method | 扫描文件或 XML 内容，返回根 namespace 和 `<entry>` 全名，为解析与修正命名空间做准备。 |

---

### `UniProtLoader.__init__(xml_input: Union[Path, str, bytes])`

**功能说明**

- 初始化解析器，确保路径有效（若输入是 `Path`），探测 namespace（调用 `detect_namespace`），并创建 `XmlParser` 以 `Entry` 为目标类型。

**参数**

- `xml_input`：支持 `Path`（必须指向文件）、XML 字符串或 XML 字节流；字符串会被传入 `BytesIO` 以便 `lxml` 读取。

**返回类型**

- `UniProtLoader` 实例。

**异常**

- `ValueError`：`Path` 输入指定位置不存在或输入字符串为空。
- `TypeError`：输入类型不是 `Path`/`str`/`bytes`。

**注意事项**

- 构造后，可反复调用 `iterate_entries` 获取新的迭代器；`self.__xml_input` 保留原始输入。

**示例**

```python
from pathlib import Path
from pybiotech.loaders.uniprot_loader import UniProtLoader

loader = UniProtLoader(Path("data/uniprot_sample.xml"))
```

---

### `UniProtLoader.iterate_entries() -> Iterator[Entry]`

**功能说明**

- 提供一个延迟生成器：如果传入的是文件，使用 `lxml.etree.iterparse` 逐条读取 `<entry>`；如果传入的是字符串，先解析整个 XML，再单独提取 entry 元素并用 `XmlParser` 建立 `Entry` 实例。

**参数**

- 无，使用构造时保存的 `xml_input`。

**返回类型**

- `Iterator[Entry]`，按文档顺序产出 `Entry` 实例；若解析失败则跳过当前节点。

**异常**

- `ValueError`：若构造阶段未能探测到 `<entry>`，`detect_namespace` 会抛出异常，从而无法正常创建加载器。

**注意事项**

- 字符串输入会在每个元素处重写 namespace，使 `xsdata` 解析器始终收到 `http://uniprot.org/uniprot`，避免前端 namespace 变化。
- 处理文件时，通过 `elem.clear()` 和删除父元素的前置节点释放内存，适合大文件。

**示例**

```python
for entry in loader.iterate_entries():
    print(entry.accession)
```

---

### `UniProtLoader.detect_namespace(xml_input: Union[Path, str, bytes]) -> Tuple[str, str]`

**功能说明**

- 扫描输入中的 `<entry>` 元素，返回实际的 namespace URI 和 `<entry>` 的完整标签；供构造函数和字符串解析流程修正 namespace。

**参数**

- `xml_input`：同构造函数，支持文件、XML 字符串或字节；字符串会尝试判断是否为路径。

**返回类型**

- `Tuple[str, str]`：`(actual_ns, entry_tag)`，其中 `entry_tag` 形式为 `"{http://...}entry"`。

**异常**

- `TypeError`：输入既不是路径也不是 XML 内容。
- `ValueError`：遍历过程中未找到任何 `<entry>` 元素。

**注意事项**

- 使用 `etree.iterparse(..., recover=True, huge_tree=True)` 实现容错；在找到目标后立即返回以缩短扫描时间。

**示例**

```python
actual_ns, entry_tag = UniProtLoader.detect_namespace("<uniprot xmlns='http://uniprot.org/uniprot'><entry/></uniprot>")
```
