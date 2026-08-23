"""
Unit tests for the ORCA input-block writer classes in
``chemsmart.io.orca.blocks``.

Each block's ``write`` method is exercised directly against an in-memory
``io.StringIO`` buffer, checking the exact text written for each
optional parameter combination.
"""

import io

from chemsmart.io.orca.blocks import (
    MaxCoreORCABlock,
    MDCIORCABlock,
    ORCABlock,
    PALORCABlock,
    SCFORCABlock,
)


class TestORCABlockBase:
    def test_stores_block_name(self):
        block = PALORCABlock(block_name="pal")
        assert block.block_name == "pal"

    def test_abstract_write_raises_not_implemented(self):
        block = ORCABlock(block_name="generic")
        try:
            block.write(io.StringIO())
        except NotImplementedError:
            pass
        else:
            raise AssertionError("Expected NotImplementedError")


class TestPALORCABlock:
    def test_write_with_nprocs(self):
        buf = io.StringIO()
        PALORCABlock(block_name="pal").write(buf, nprocs=8)
        assert buf.getvalue() == "%pal nprocs 8 end\n\n"

    def test_write_without_nprocs_is_noop(self):
        buf = io.StringIO()
        PALORCABlock(block_name="pal").write(buf)
        assert buf.getvalue() == ""


class TestMaxCoreORCABlock:
    def test_write_with_maxcore(self):
        buf = io.StringIO()
        MaxCoreORCABlock(block_name="maxcore").write(buf, maxcore=2000)
        assert buf.getvalue() == "%maxcore 2000\n\n"

    def test_write_without_maxcore_is_noop(self):
        buf = io.StringIO()
        MaxCoreORCABlock(block_name="maxcore").write(buf)
        assert buf.getvalue() == ""


class TestSCFORCABlock:
    def test_write_convergence(self):
        buf = io.StringIO()
        SCFORCABlock(block_name="scf").write(buf, convergence="tight")
        assert "convergence tight" in buf.getvalue()

    def test_write_rotate(self):
        buf = io.StringIO()
        SCFORCABlock(block_name="scf").write(buf, rotate="{48, 49, 90, 1, 1}")
        assert "rotate {48, 49, 90, 1, 1} end" in buf.getvalue()

    def test_write_dryrun(self):
        buf = io.StringIO()
        SCFORCABlock(block_name="scf").write(buf, dryrun=True)
        assert "DryRun true" in buf.getvalue()

    def test_write_scfmeminfo(self):
        buf = io.StringIO()
        SCFORCABlock(block_name="scf").write(buf, scfmeminfo=True)
        assert "Print[P_SCFMemInfo] 1" in buf.getvalue()

    def test_write_guessmode(self):
        buf = io.StringIO()
        SCFORCABlock(block_name="scf").write(buf, guessmode="CMatrix")
        assert "GuessMode CMatrix" in buf.getvalue()

    def test_write_autotrah(self):
        buf = io.StringIO()
        SCFORCABlock(block_name="scf").write(buf, autotrah=True)
        assert "AutoTRAH true" in buf.getvalue()

    def test_write_flipspin_and_finalms_together(self):
        buf = io.StringIO()
        SCFORCABlock(block_name="scf").write(
            buf, flipspin="17,38,56", finalms=0.0
        )
        assert "FlipSpin 17,38,56 FinalMs 0.0" in buf.getvalue()

    def test_flipspin_without_finalms_is_ignored(self):
        buf = io.StringIO()
        SCFORCABlock(block_name="scf").write(buf, flipspin="17,38,56")
        assert "FlipSpin" not in buf.getvalue()

    def test_write_defaults_is_noop(self):
        buf = io.StringIO()
        SCFORCABlock(block_name="scf").write(buf)
        assert buf.getvalue() == ""


class TestMDCIORCABlock:
    def test_write_minimal_defaults(self):
        buf = io.StringIO()
        MDCIORCABlock(block_name="mdci").write(buf)
        content = buf.getvalue()
        assert content.startswith("%mdci\n")
        assert content.endswith("end\n\n")
        # default-valued options should still appear
        assert "pCCSDAB 0.1" in content
        assert "STol 1e-05" in content
        assert "LShift 0.3" in content
        assert "MaxDIIS 7" in content

    def test_write_all_options(self):
        buf = io.StringIO()
        MDCIORCABlock(block_name="mdci").write(
            buf,
            citype="CCSD",
            ewin=(-3, 1000),
            singles=True,
            triples=1,
            brueckner=True,
            denmat="orbopt",
            zsimple=True,
            useqros=True,
            localize="PM",
            natorbiters=0,
            pccsdab=0.5,
            pccsdcd=0.6,
            pccsdef=0.7,
            kcopt="KC_RI",
            printLevel=3,
            maxiter=50,
            maxcore=1000,
            trafotype="trafo_ri",
            stol=1e-6,
            lshift=0.2,
            maxdiis=5,
            incore=2,
        )
        content = buf.getvalue()
        assert "citype CCSD" in content
        assert "ewin -3,1000" in content
        assert "Singles true" in content
        assert "Triples 1" in content
        assert "Brueckner true" in content
        assert "Denmat orbopt" in content
        assert "ZSimple true" in content
        assert "UseQROs" in content
        assert "Localize PM" in content
        assert "NatOrbIters 0" in content
        assert "pCCSDAB 0.5" in content
        assert "pCCSDCD 0.6" in content
        assert "pCCSDEF 0.7" in content
        assert "KCOpt KC_RI" in content
        assert "PrintLevel 3" in content
        assert "MaxIter 50" in content
        assert "MaxCore 1000" in content
        assert "TrafoType trafo_ri" in content
        assert "STol 1e-06" in content
        assert "LShift 0.2" in content
        assert "MaxDIIS 5" in content
        assert "InCore 2" in content

    def test_write_falsy_flags_omitted(self):
        buf = io.StringIO()
        MDCIORCABlock(block_name="mdci").write(
            buf,
            singles=False,
            brueckner=False,
            zsimple=False,
            useqros=False,
        )
        content = buf.getvalue()
        assert "Singles" not in content
        assert "Brueckner" not in content
        assert "ZSimple" not in content
        assert "UseQROs" not in content
