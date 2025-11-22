# SDF Loader

模块 `pybiotech.loaders.sdf_loader` 提供一组针对 SDF 文件/目录/文本的读入工具，围绕 `RDKit` 的 `Mol` 对象完成批量筛选、错误忽略和日志控制。

## 模块概览

| 名称 | 类型 | 描述 |
| --- | --- | --- |
| `SDFLoader` | class | 聚合所有加载逻辑的载体类，以静态方法暴露行为，无需实例化。 |
| `SDFLoader.closeWarning` | static method | 关闭 RDKit 的日志输出，仅保留错误级别。 |
| `SDFLoader.openWarning` | static method | 恢复 RDKit 的默认（WARNING）日志级别。 |
| `SDFLoader.readDataFromFile` | static method | 按索引范围从单个 `.sdf` 文件读取分子，返回有效分子和计数。 |
| `SDFLoader.readDataFromDir` | static method | 遍历目录（可递归）中所有 `.sdf` 文件，批量读取分子。 |
| `SDFLoader.readDataFromText` | static method | 通过临时文件把 SDF 文本内容转为文件后再读取，复用 `readDataFromFile`。 |

---

### `SDFLoader.closeWarning()`

**功能说明**

- 直接把 `rdkit.RDLogger` 的等级调为 `ERROR`，用于在批量读取或批处理脚本里抑制冗长的 RDKit 警告输出。

**参数**

- 无。

**返回类型**

- `None`

**异常**

- 无。

**注意事项**

- 本方法在同一进程中的其他 RDKit 调用也会感知到日志等级变化；需要 `openWarning` 恢复。

**示例**

```python
from pybiotech.loaders.sdf_loader import SDFLoader

SDFLoader.closeWarning()
```

---

### `SDFLoader.openWarning()`

**功能说明**

- 把 `RDLogger` 的等级恢复为 `WARNING`，确保后续需要查看警告时可见。

**参数**

- 无。

**返回类型**

- `None`

**异常**

- 无。

**注意事项**

- 配对使用 `closeWarning` 时，建议在长时间无须追踪警告后再调用；对短脚本而言不用关闭即可。

**示例**

```python
SDFLoader.openWarning()
```

---

### `SDFLoader.readDataFromFile(sdfDoc: str, startIndex: int = 0, endIndex: int = -1, ignore_error: bool = True)`

**功能说明**

- 读取单个 SDF 文件中 `[startIndex, endIndex)` 范围内的结构，跳过空分子；返回有效分子列表、理论应读数量与实际读取数量，便于后续 QA。

**参数**

- `sdfDoc`：必须存在的文件路径；使用 `Chem.ForwardSDMolSupplier` 以 `sanitize=False` 读取原子。
- `startIndex`：起始位置（包含），默认 0，不能为负数。
- `endIndex`：结束位置（不包含），默认 `-1` 表示读取到文件尾。
- `ignore_error`：若为 `True`，遇到 `None` 或无原子的记录会静默跳过；否则抛 `ValueError`。

**返回类型**

- `tuple[list[Mol], int, int]`：有效分子列表；期望读取数量（`endIndex-startIndex` 或剩余到末尾）；实际有效读取数量。

**异常**

- `FileNotFoundError`：`sdfDoc` 不存在。
- `ValueError`：`startIndex < 0`、`endIndex != -1 且 endIndex <= startIndex`、遇到无效分子且 `ignore_error=False`。

**注意事项**

- 使用 `removeHs=False` 保留显式氢以便后续处理；`sanitize=False` 可避免提前抛错。若需要 `Chem.AddHs` 可在后续 pipeline 完成。
- `endIndex` 设置为 `-1` 时，函数会遍历整个文件，并计算 `totalCount` 以返回 `toReadCount`。

**示例**

```python
from rdkit.Chem import Mol

mols, should_read, actual = SDFLoader.readDataFromFile(
    "data/biologics.sdf",
    startIndex=10,
    endIndex=50,
    ignore_error=False,
)
print(f"期望{should_read}条，已收集{actual}条")
```

---

### `SDFLoader.readDataFromDir(sdfDir: str, recursive: bool = True, ignore_error: bool = True)`

**功能说明**

- 遍历给定目录下所有 `.sdf` 文件（可递归），逐个调用 `readDataFromFile` 合并所有结构。

**参数**

- `sdfDir`：必须存在且为目录。
- `recursive`：是否递归子目录，默认为 `True`。
- `ignore_error`：传给 `readDataFromFile` 的 `ignore_error`，决定是否在遇到异常时跳过。

**返回类型**

- `tuple[list[Mol], int, int]`：包含所有文件的有效分子列表、按索引理应读取的总数和实际有效分子数。

**异常**

- `FileNotFoundError`：目录不存在。
- `TypeError`：传入的 `sdfDir` 不是目录。

**注意事项**

- `glob.glob(..., recursive=recursive)` 用于匹配任意深度的 `.sdf`，`recursive=False` 会限制在当前目录。
- 如果目录下无 `.sdf` 文件，返回的列表为空、计数为 0。

**示例**

```python
mols, expected, collected = SDFLoader.readDataFromDir(
    "data/sdf",
    recursive=True,
    ignore_error=True,
)
```

---

### `SDFLoader.readDataFromText(sdfText: str, ignore_error: bool = True)`

**功能说明**

- 通过临时文件将输入文本写入磁盘，复用 `readDataFromFile` 的逻辑来解析 SDF 格式。

**参数**

- `sdfText`：包含完整 SDF 内容的文本字符串。
- `ignore_error`：是否忽略无效记录，默认 `True`。

**返回类型**

- `tuple[list[Mol], int, int]`：同 `readDataFromFile`。

**异常**

- 异常直接透传自 `readDataFromFile`，如文件缺失或索引错误。

**注意事项**

- 临时文件会写入磁盘并在读取后尝试删除；若删除失败会打印异常。
- 避免传入空字符串或非 SDF 内容否则 `Chem.ForwardSDMolSupplier` 可能返回全 `None`。

**示例**

```python
sdf_content = open("data/example.sdf").read()
mols, should_read, actual = SDFLoader.readDataFromText(sdf_content)
```
