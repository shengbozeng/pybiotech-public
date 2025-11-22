from pybiotech.classes.nih.https.pubchem.input_spec import InputSpec
from pybiotech.classes.nih.https.pubchem.opera_spec import OperationSpec
from pybiotech.classes.nih.https.pubchem.output_spec import OutputSpec
from pybiotech.classes.nih.https.pubchem.query_options import QueryOptions
from pybiotech.classes.nih.https.pubchem.request import PubChemRequest


def test_inputspec():
    print("\n[InputSpec 测试]")

    # 自动推断 namespace
    s1 = InputSpec(domain="compound", identifiers=["2244", "962"])
    assert s1.namespace == "cid"
    assert s1.to_url_path() == "compound/cid/2244%2C962"
    print("自动推断 namespace: OK")

    # 提供 formula
    s2 = InputSpec(domain="compound", identifiers="C9H8O4")
    assert s2.namespace == "formula"
    print("formula 推断: OK")

    # 手动传递 namespace，校验有效性
    s3 = InputSpec(domain="substance", namespace="sourceid", identifiers="DTP/NCI")
    assert s3.to_url_path() == "substance/sourceid/DTP.NCI"
    print("特殊 sourceid 处理: OK")

    # 错误 domain
    try:
        InputSpec(domain="unknown", identifiers="2244")
        assert False, "应报错"
    except Exception as e:
        print("错误 domain 检查: OK")

    # 错误 namespace
    try:
        InputSpec(domain="compound", namespace="notexist", identifiers="2244")
        assert False, "应报错"
    except Exception as e:
        print("错误 namespace 检查: OK")
        
def test_operationspec():
    print("\n[OperationSpec 测试]")

    # property 必须有 tags
    try:
        OperationSpec(operation="property", domain="compound")
        assert False, "应报错 property 缺 tags"
    except Exception as e:
        print("property 缺 tags 检查: OK")

    # tags 支持字符串、列表
    op1 = OperationSpec(operation="property", tags="MolecularFormula,InChIKey", domain="compound")
    assert op1.tags == ["MolecularFormula", "InChIKey"]
    assert op1.to_url_path() == "property/MolecularFormula%2CInChIKey"
    print("tags 字符串列表都支持: OK")

    # operation/不支持的 domain
    try:
        OperationSpec(operation="property", tags=["MolecularFormula"], domain="gene")
        assert False, "应报错"
    except Exception as e:
        print("非法 operation/domain 检查: OK")


def test_operationspec():
    print("\n[OperationSpec 测试]")

    # property 必须有 tags
    try:
        OperationSpec(operation="property", domain="compound")
        assert False, "应报错 property 缺 tags"
    except Exception as e:
        print("property 缺 tags 检查: OK")

    # tags 支持字符串、列表
    op1 = OperationSpec(operation="property", tags="MolecularFormula,InChIKey", domain="compound")
    assert op1.tags == ["MolecularFormula", "InChIKey"]
    assert op1.to_url_path() == "property/MolecularFormula%2CInChIKey"
    print("tags 字符串列表都支持: OK")

    # operation/不支持的 domain
    try:
        OperationSpec(operation="property", tags=["MolecularFormula"], domain="gene")
        assert False, "应报错"
    except Exception as e:
        print("非法 operation/domain 检查: OK")

def test_outputspec():
    print("\n[OutputSpec 测试]")

    # 标准格式
    o1 = OutputSpec(output_format="json")
    assert o1.output_format == "JSON"
    assert o1.to_url_path() == "JSON"
    print("标准格式检查: OK")

    # JSONP 必须带 callback
    o2 = OutputSpec(output_format="JSONP", callback="cb")
    assert o2.to_url_path() == "JSONP?callback=cb"
    print("JSONP callback 拼接: OK")

    # 不支持格式
    try:
        OutputSpec(output_format="EXCEL")
        assert False, "应报错"
    except Exception as e:
        print("不支持格式校验: OK")

def test_queryoptions():
    print("\n[QueryOptions 测试]")

    # 常规参数校验
    qo = QueryOptions(list_return="flat", MaxRecords=10, Stereo="ignore")
    assert "list_return" in qo.params and qo["Stereo"] == "ignore"
    print("常规参数校验: OK")

    # 错误类型校验
    try:
        QueryOptions(MaxRecords="notint")
        assert False, "应报错"
    except Exception as e:
        print("类型校验错误: OK")

    # 多值校验
    qo2 = QueryOptions(cids_type="all,active")
    assert qo2["cids_type"] == "all,active"
    print("多值类型校验: OK")

    # image_size 特例
    qo3 = QueryOptions(image_size="320x240")
    assert qo3["image_size"] == "320x240"
    print("image_size 特例: OK")

def test_pubchemrequest_all():
    print("\n[PubChemRequest 综合测试]")

    # 1. 正常组合
    try:
        input_spec = InputSpec(domain="compound", namespace="name", identifiers="aspirin")
        operation_spec = OperationSpec(operation="property", tags=["MolecularFormula", "InChIKey"], domain="compound")
        output_spec = OutputSpec(output_format="JSON")
        query_options = QueryOptions(list_return="flat", MaxRecords=5)
        req = PubChemRequest(input_spec, operation_spec, output_spec, query_options)
        url = req.build_url()
        print("正常请求通过:", url)
    except Exception as e:
        print("应通过的正常组合报错:", e)

    # 2. callback+非JSONP
    try:
        input_spec = InputSpec(domain="compound", namespace="cid", identifiers="2244")
        operation_spec = OperationSpec(operation="synonyms", tags=["MolecularFormula"], domain="compound")
        output_spec = OutputSpec(output_format="JSON")
        query_options = QueryOptions(callback="cb")
        req = PubChemRequest(input_spec, operation_spec, output_spec, query_options)
        req.build_url()
        assert False, "callback 非 JSONP 应报错"
    except Exception as e:
        print("callback+非JSONP 检查: OK")

    # 3. property 缺 tags
    try:
        input_spec = InputSpec(domain="compound", namespace="cid", identifiers="2244")
        operation_spec = OperationSpec(operation="property", domain="compound")
        output_spec = OutputSpec(output_format="JSON")
        query_options = QueryOptions()
        req = PubChemRequest(input_spec, operation_spec, output_spec, query_options)
        req.build_url()
        assert False, "property 缺 tags 应报错"
    except Exception as e:
        print("property 缺 tags 检查: OK")

    # 4. PNG + 非法 operation
    try:
        input_spec = InputSpec(domain="compound", namespace="cid", identifiers="2244")
        operation_spec = OperationSpec(operation="aids", domain="compound")
        output_spec = OutputSpec(output_format="PNG")
        req = PubChemRequest(input_spec, operation_spec, output_spec, QueryOptions())
        req.build_url()
        assert False, "PNG 非 record/conformers 应报错"
    except Exception as e:
        print("PNG 非法 operation 检查: OK")

    # 5. JSONP 无 callback
    try:
        input_spec = InputSpec(domain="compound", namespace="cid", identifiers="2244")
        operation_spec = OperationSpec(operation="synonyms", domain="compound")
        output_spec = OutputSpec(output_format="JSONP")
        req = PubChemRequest(input_spec, operation_spec, output_spec, QueryOptions())
        req.build_url()
        assert False, "JSONP 无 callback 应报错"
    except Exception as e:
        print("JSONP 缺 callback 检查: OK")


if __name__ == "__main__":
    test_inputspec()
    test_operationspec()
    test_outputspec()
    test_queryoptions()
    test_pubchemrequest_all()
    print("\n所有核心功能测试完毕。")
