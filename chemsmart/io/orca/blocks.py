import logging
from abc import abstractmethod

logger = logging.getLogger(__name__)


class ORCABlock:
    """
    Abstract base class for ORCA input file blocks.
    
    This class provides the interface for handling different types of ORCA
    input blocks. Each block type should inherit from this class and implement
    the write method to output the block content to an input file.
    """
    
    def __init__(self, block_name):
        """
        Initialize ORCA block object.
        
        Args:
            block_name (str): Name of the ORCA input block to handle
        """
        self.block_name = block_name

    @abstractmethod
    def write(self, f, **kwargs):
        """
        Abstract method to write block content to file.
        
        Args:
            f: File object to write to
            **kwargs: Block-specific parameters
            
        Raises:
            NotImplementedError: Must be implemented by subclasses
        """
        raise NotImplementedError


class PALORCABlock(ORCABlock):
    """
    ORCA input %pal block for parallel processing specification.
    
    This block controls the number of processors used for ORCA calculations.
    Essential for parallel computation setup in quantum chemistry calculations.
    """

    def __init__(self, block_name):
        """
        Initialize PAL ORCA block.
        """
        super().__init__(block_name=block_name)

    def write(
        self,
        f,
        nprocs=None,  # pal
    ):
        """
        Write %pal block to ORCA input file.
        
        Args:
            f: File object to write to
            nprocs (int, optional): Number of processors to use for calculation
        """
        if nprocs is not None:
            f.write(f"%pal nprocs {nprocs} end\n\n")


class MaxCoreORCABlock(ORCABlock):
    """
    ORCA input %maxcore block for maximum memory per core specification.
    
    This block controls the maximum amount of memory allocated per core
    for ORCA calculations. Important for managing memory usage in large
    quantum chemistry calculations.
    """

    def __init__(self, block_name):
        """
        Initialize MaxCore ORCA block.
        
        Args:
            block_name (str): Name of the block (typically 'maxcore')
        """
        super().__init__(block_name=block_name)

    def write(
        self,
        f,
        maxcore=None,
    ):
        """
        Write %maxcore block to ORCA input file.
        
        Args:
            f: File object to write to
            maxcore (int, optional): Maximum memory per core in MB
        """
        if maxcore is not None:
            f.write(f"%maxcore {maxcore}\n\n")


class SCFORCABlock(ORCABlock):
    """
    ORCA input %scf block for SCF algorithm parameters and options.
    
    This block controls various Self-Consistent Field (SCF) calculation
    parameters including convergence criteria, orbital manipulation,
    and broken-symmetry solutions for open-shell systems.
    """

    def __init__(self, block_name):
        """
        Initialize SCF ORCA block.
        
        Args:
            block_name (str): Name of the block (typically 'scf')
        """
        super().__init__(block_name=block_name)

    def write(
        self,
        f,
        convergence=None,  # convergence tight <- requires `convergence='tight'`
        rotate=None,  # rotate {48, 49, 90, 1, 1} end <- requires `rotate='{48, 49, 90, 1, 1}'`
        dryrun=False,  # dryrun true <- requires `dryrun=True`
        scfmeminfo=False,  # Print[P_SCFMemInfo] 1 <- requires `scfmeminfo=True`
        guessmode=None,  # GuessMode CMatrix  <- requires `guessmode='CMatrix'`
        autotrah=False,  # AutoTRAH true
        # for generating broken-symmetry solutions:
        flipspin=None,  # FlipSpin 1
        # Flip spin is a vector and you can give a list of atom on which you want to have the spin
        # flipped.  For example
        # FlipSpin 17,38,56
        # REMEMBER: counting starts at atom 0!
        finalms=None,  # FinalMs 0.0
        # The desired Ms value of the broken symmetry determinant.
        # This value MUST be given since the program cannot determine it by itself
    ):
        """
        Write %scf block to ORCA input file.
        """
        if convergence is not None:
            f.write(f"%scf\n  convergence {convergence} \nend\n\n")
        if rotate is not None:
            f.write(f"%scf\n  rotate {rotate} end \nend\n\n")
        if dryrun:
            f.write("%scf\n  DryRun true \nend\n\n")
        if scfmeminfo:
            f.write("%scf\n   Print[P_SCFMemInfo] 1\n end\n\n")
        if guessmode is not None:
            f.write(f"%scf\n  GuessMode {guessmode}\n end\n\n")
        if autotrah:
            f.write("%scf\n  AutoTRAH true \nend\n\n")
        if (
            flipspin is not None and finalms is not None
        ):  # flipSpin and finalMs values need to be specified together
            f.write(f"%scf\n  FlipSpin {flipspin} FinalMs {finalms}\n end\n\n")


class MDCIORCABlock(ORCABlock):
    """
    ORCA input %mdci block for Matrix Driven CI methods.
    
    This block controls various correlation methods including CCSD, CISD,
    QCISD, and CEPA methods. It provides extensive options for controlling
    convergence, memory usage, and calculation details.
    """

    def __init__(self, block_name):
        """
        Initialize MDCI ORCA block.
        
        Args:
            block_name (str): Name of the block (typically 'mdci')
        """
        super().__init__(block_name=block_name)

    def write(  # noqa: PLR0912
        self,
        f,
        citype=None,
        ewin=None,  # need a 2-tuple
        singles=False,
        triples=None,  # options: 0, 1, 2
        brueckner=False,
        denmat=None,  # options: none, linearized, unrelaxed, orbopt
        zsimple=False,
        useqros=False,
        localize=None,  # options: 0, PM, FB
        natorbiters=None,
        pccsdab=0.1,
        pccsdcd=0.1,
        pccsdef=0.1,
        kcopt=None,  # options: KC_MO, KC_AOBLAS, KC_AO, KC_RI, KC_RI2, KC_AOX
        printLevel=None,  # options: 2, 3, 4
        maxiter=None,  # integer number
        maxcore=None,  # integer number
        trafotype=None,  # trafo_jk, trafo_ri, trafo_full
        stol=1e-5,
        lshift=0.3,
        maxdiis=7,
        incore=None,  # options: 0,1,2,3,4,5
    ):
        """Write various possible options in MDCI module in ORCA.

        Arguments:
            f:      file object to write to
            citype: options:
                        CISD            # CI singles+doubles
                        QCISD           # quadratic CI (singles+doubles)
                        CCSD            # coupled-cluster singles+doubles
                        CEPA_1          # coupled-electron pair approximation '1'
                        CEPA_2          # coupled-electron pair approximation '2'
                        CEPA_3          # coupled-electron pair approximation '3'
                        CPF_1           # Coupled-pair functional approximation '1'
                        CPF_2           # (note that CPF/1 is identical with
                        CPF_3           #  the original CPF of Ahlrichs et al.)
                        VCEPA_1         # Variational CEPA approximation '1'
                        VCEPA_2         #
                        VCEPA_3         #
                        NCEPA_1         # our slightly modified versions of CEPA
                        NCEPA_2         # and CPF
                        NCEPA_3         #
                        NCPF_1          #
                        NCPF_2          #
                        NCPF_3          #
                        VNCEPA_1        #
                        VNCEPA_2        #
                        VNCEPA_3        #
                        ACPF            # averaged coupled-pair functional
                        ACPF_2          # Gdanitz modification of it
                        NACPF           # our modification of it
                        AQCC            # Szalay + Bartlett
                        SEOI            # a strictly size extensive energy functional maintaining unitary invariance
                                          (not yet available for UHF)
                        MP3             # MP3 calculation. With UseSCS=true it is SCS-MP3
            ewin: -3,1e3                # orbital energy window to determine which MOs are included in the treatment
            singles: options:
                        true            # include single excitations in the treatment (default true)
                        false
            triples: options:
                        0               # (T) correction in CCSD(T)/QCISD(T) default is no triples
                        1               # Algorithm 1 (lots of memory, fast)
                        2               # Algorithm 2 (less memory, about 2x slower)
            brueckner: options:
                        true            # use Brueckner orbitals (default false)
                        false
            denmat: options:
                        none            # no evaluation of density matrix
                        linearized      # density matrix obtained by retaining only CEPA_0-like terms, i.e.,
                                        # those linear in the excitation amplitudes
                        unrelaxed       # unrelaxed density matrices, i.e., density matrices without orbital relaxation
                        orbopt          # perform orbital optimization yielding fully relaxed density matrices
                                        # (if citype chosen as CCSD or QCISD this option implies evaluation of the Z
                                          vector).
                                        # (default: linearized)
            zsimple: options:
                        true            # simplified evaluation of the Z vector in case of orbital optimized CCD
                                        # (citype chosen as CCSD or QCISD and Denmat as orbopt) by using an analytical
                                           formula
                        false           # explicit solution of Z vector equations in case of orbital optimized CCD
                                        # (default: false)
            useqros:                     # use of quasi-restricted orbitals (default false)
            localize: options:
                        0               # use localized MOs. Presently very little use is made of locality. It may help
                                        # for interpretations. Localization is incompatible with the (T) correction
                        PM              # Use Pipek-Mezey localized MOs
                        FB              # Use Foster-Boys localized MOs
            natorbiters: options:
                        0               # Perform natural orbital iterations. default is none. Not possible for CCSD
                                          and QCISD
            pccsdab:                     # the three parameters for parametrized
            pccsdcd:                     # coupled-cluster (default is 1.0 which corresponds to normal CCSD)
            pccsdef:                     # this defines how the rate limiting step is handled
                                        # MO and AOX need lots of disk and I/O but if they can be done they are fast
            kcopt: options:
                    KC_MO               # Perform full 4-index transformation
                    KC_AOBLAS           # AO direct with BLAS (preferred)
                    KC_AO               # AO direct handling of 3,4 externals
                    KC_RI               # RI approximation of 3,4 externals
                    KC_RI2              # Alternative RI (not recommended)
                    KC_AOX              # Do it from stored AO exchange integrals
            printLevel: options:
                    2                   # Control the amount of output.
                    3                   # For 3 and higher things like pair correlation energies are printed.
            maxiter: 35                  # Max. number of iterations
            trafotype: options:         # How the integral transformation is done.
                                        # Note that it is fine to do AOX or AO or AOBLAS together with trafo_ri
                    trafo_jk            # Partial trafo to J+K operators
                    trafo_ri            # RI transformation of all integrals up to 2-externals (3-ext for (T)) and
                                          rest on the fly
                    trafo_full          # Full four index transformation.  Automatically chosen for KCOpt=KC_MO
            maxcore: 350                # Memory in MB - used for integral trafos and batching and for storage of
                                        # integrals and amplitudes # don't be too generous
            stol: 1e-5                   # Max. element of the residual vector for convergence check
            lshift: 0.3                  # Level shift to be used in update of coefficients
            maxdiis:  7                  # Max number of DIIS vectors to be stored.
                                        # this lets you control how much and what is residing in central memory.
                                        # May speed up things. Note that MaxCore is not respected here 9
            incore: options:
                    0                   # nothing in core
                    1                   #  + sigma-vector and amplitudes (default)
                    2                   #  + Jij(a,b) Kij(a,b) operators
                    3                   #  + DIIS vectors
                    4                   #  + 3-exernal integral Kia(b,c)
                    5                   #  + 4-external integrals Kab(c,d)
                                        # this is identical to ALL the default is AUTO which means that incore is
                                          chosen based on MaxCore.
        """
        f.write("%mdci\n")

        if citype:
            f.write(f"citype {citype}\n")
        if ewin:
            f.write(
                f"ewin {ewin[0]},{ewin[1]}\n # orbital energy window to determine which\n"
                f"                                     # MOs are included in the treatment\n"
            )
        if singles:
            f.write("Singles true\n")
        if triples is not None:
            f.write(f"Triples {triples}\n")
        if brueckner:
            f.write("Brueckner true\n")
        if denmat is not None:
            f.write(f"Denmat {denmat}\n")
        if zsimple:
            f.write("ZSimple true\n")
        if useqros:
            f.write("UseQROs\n")
        if localize is not None:
            f.write(f"Localize {localize}\n")
        if natorbiters is not None:
            f.write(f"NatOrbIters {natorbiters}\n")
        if pccsdab is not None:
            f.write(f"pCCSDAB {pccsdab}\n")
        if pccsdcd is not None:
            f.write(f"pCCSDCD {pccsdcd}\n")
        if pccsdef is not None:
            f.write(f"pCCSDEF {pccsdef}\n")
        if kcopt is not None:
            f.write(f"KCOpt {kcopt}\n")
        if printLevel is not None:
            f.write(f"PrintLevel {printLevel}\n")
        if maxiter is not None:
            f.write(f"MaxIter {maxiter}\n")
        if maxcore is not None:
            f.write(f"MaxCore {maxcore}\n")
        if trafotype is not None:
            f.write(f"TrafoType {trafotype}\n")
        if stol is not None:
            f.write(f"STol {stol}\n")
        if lshift is not None:
            f.write(f"LShift {lshift}\n")
        if maxdiis is not None:
            f.write(f"MaxDIIS {maxdiis}\n")
        if incore is not None:
            f.write(f"InCore {incore}\n")
            
        f.write("end\n\n")
