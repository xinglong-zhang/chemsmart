# Excluded engine artefacts

This tree carries the complete human-readable record of both campaigns. It
deliberately omits ORCA's binary restart artefacts, which are large,
unreadable, regenerable from the retained inputs, and read by nothing in
`REVIEW.md`, `PEER-REVIEW.md`, or the typed analysis plane:

| extension | what it is | why omitted |
|---|---|---|
| `.gbw` | converged molecular-orbital coefficients | scales as (basis functions)^2; 680-function def2-TZVP saddles give ~7 MB each |
| `.densities`, `.densitiesinfo` | electron-density containers | derived from the `.gbw` |
| `.opt` | binary optimisation-restart state | supersedes nothing in the `.out` trajectory |
| `.cpcm`, `.cpcm_corr` | continuum-solvation cavity state | the applied solvation terms are printed in the `.out` |

Every file below is listed with the SHA-256 of the bytes that were on disk
at publication, so a reader who obtains one can verify it belongs to this
run, and a claim that depends on one can be traced to a file that is named
rather than absent. What is retained is sufficient for every number in the
reviews: the engine logs (`.out`), the exact submitted inputs (`.inp`), the
Hessians (`.hess`), the ORCA property dumps (`_property.txt`), all
geometries (`.xyz`), and the agent's own receipts, event streams, project
YAML and retained preview inputs.

| sha256 | bytes | path |
|---|---|---|
| `5bb5b513840e969655c4fb0e68c93a3deb1a03f2707633e37b72853c3f97a2f2` | 537715 | `ino3-r12/workspace/nodes/cat-sp-svp-b3lyp/geom-pme3-cat-dbl-min_sp_sp_cpcm_acetonitrile.cpcm` |
| `34c605d988a0c8be6897b803e040ee79c55fb677b22c049b6af1b8096bb030bb` | 54132 | `ino3-r12/workspace/nodes/cat-sp-svp-b3lyp/geom-pme3-cat-dbl-min_sp_sp_cpcm_acetonitrile.cpcm_corr` |
| `0532ed5021af34552c554c0baec6a005b8a99d345a717f5d66c5c8cc35d181fc` | 3591200 | `ino3-r12/workspace/nodes/cat-sp-svp-b3lyp/geom-pme3-cat-dbl-min_sp_sp_cpcm_acetonitrile.densities` |
| `9c3e726dfb6ceb90f7200e26cb7c2008307eb39e46b6ab4cb8ff6be58715538d` | 3160 | `ino3-r12/workspace/nodes/cat-sp-svp-b3lyp/geom-pme3-cat-dbl-min_sp_sp_cpcm_acetonitrile.densitiesinfo` |
| `34ea14a876925918203e5ff2cc28ef668939d1396ccb7cdeb6d9f42cf8cf7409` | 3107508 | `ino3-r12/workspace/nodes/cat-sp-svp-b3lyp/geom-pme3-cat-dbl-min_sp_sp_cpcm_acetonitrile.gbw` |
| `84ef553f9f2ee6561decb88a817e8e3d4faffeff50c7cd7fb06197f138ac9bf2` | 537715 | `ino3-r12/workspace/nodes/cat-sp-svp-pbe0/geom-pme3-cat-dbl-min_sp_sp_cpcm_acetonitrile.cpcm` |
| `4f03055fbfc8e1c56063c947b2dd71d26ae471fe85506c76ca4e5eea86ac53be` | 54132 | `ino3-r12/workspace/nodes/cat-sp-svp-pbe0/geom-pme3-cat-dbl-min_sp_sp_cpcm_acetonitrile.cpcm_corr` |
| `acb8ef0978a3eb8c822c058268edd844b049ab8670c29a25f349d9eeab1cf6b3` | 3591200 | `ino3-r12/workspace/nodes/cat-sp-svp-pbe0/geom-pme3-cat-dbl-min_sp_sp_cpcm_acetonitrile.densities` |
| `5ba087b55fd3a3d6d781d5a9f455c266193dee3db4f52da84658a13aadc05d46` | 3160 | `ino3-r12/workspace/nodes/cat-sp-svp-pbe0/geom-pme3-cat-dbl-min_sp_sp_cpcm_acetonitrile.densitiesinfo` |
| `92d2c0a601ebe7c293b6452890ae2de1bbcb6083427ef167ebea876551195265` | 3107508 | `ino3-r12/workspace/nodes/cat-sp-svp-pbe0/geom-pme3-cat-dbl-min_sp_sp_cpcm_acetonitrile.gbw` |
| `9cf68def9ee4bf7f895d68b8272f46f2706d0354552d80f487b76759b89be463` | 537715 | `ino3-r12/workspace/nodes/cat-sp-svp-tpssh/geom-pme3-cat-dbl-min_sp_sp_cpcm_acetonitrile.cpcm` |
| `ca2417ea31795cd81f262c01d1f3d952fd6e4e7ab719dcfdd98b94db8192fae2` | 54132 | `ino3-r12/workspace/nodes/cat-sp-svp-tpssh/geom-pme3-cat-dbl-min_sp_sp_cpcm_acetonitrile.cpcm_corr` |
| `b630e154dbccca4952996cc7d77d79a8b3c8d9e3188c17bcd03d894e30ceb46b` | 3591200 | `ino3-r12/workspace/nodes/cat-sp-svp-tpssh/geom-pme3-cat-dbl-min_sp_sp_cpcm_acetonitrile.densities` |
| `8f7320b0b420025d214d44d87ef8e71b4f6826277542409299761590eb081d77` | 3160 | `ino3-r12/workspace/nodes/cat-sp-svp-tpssh/geom-pme3-cat-dbl-min_sp_sp_cpcm_acetonitrile.densitiesinfo` |
| `c540d9ccb2d7a5eebe68ed6eb81e881854c3b405fb4ccef6000a69d2ac4076d9` | 3107508 | `ino3-r12/workspace/nodes/cat-sp-svp-tpssh/geom-pme3-cat-dbl-min_sp_sp_cpcm_acetonitrile.gbw` |
| `d7d21ea29f4e6551068ea4558d2cec84a7dfc5cffdba7ef480907487dba735a6` | 519727 | `ino3-r12/workspace/nodes/n0-sp-svp-b3lyp/geom-pme3-neutral-min_sp_sp_cpcm_acetonitrile.cpcm` |
| `4387ebf1b60f324ff599061a6ae0007e96d925f8b48aed64a9eb5b2f81be2f87` | 52313 | `ino3-r12/workspace/nodes/n0-sp-svp-b3lyp/geom-pme3-neutral-min_sp_sp_cpcm_acetonitrile.cpcm_corr` |
| `012a8d50d6fd43ccba090f347fb45d871d07c06677b09c4879acbd864db5ff39` | 1795600 | `ino3-r12/workspace/nodes/n0-sp-svp-b3lyp/geom-pme3-neutral-min_sp_sp_cpcm_acetonitrile.densities` |
| `d6b4bf82165c05fabaef6a2761ddfe4ed9014a03351f9299cdc9bdfb9f4e460b` | 1838 | `ino3-r12/workspace/nodes/n0-sp-svp-b3lyp/geom-pme3-neutral-min_sp_sp_cpcm_acetonitrile.densitiesinfo` |
| `c4ad1d7071fdc465d2cef2b3e53ad02f9d8813d3598bda1b0042c1d64794e0ed` | 2201668 | `ino3-r12/workspace/nodes/n0-sp-svp-b3lyp/geom-pme3-neutral-min_sp_sp_cpcm_acetonitrile.gbw` |
| `108089b3809be711f6290422118fa5cc6b9be6bb4f8e44032e833f50ae5d9f42` | 519727 | `ino3-r12/workspace/nodes/n0-sp-svp-pbe0/geom-pme3-neutral-min_sp_sp_cpcm_acetonitrile.cpcm` |
| `395749a1f2c22089ed80c4db08b281d21bdde549666e189f890eceb0e04fc109` | 52313 | `ino3-r12/workspace/nodes/n0-sp-svp-pbe0/geom-pme3-neutral-min_sp_sp_cpcm_acetonitrile.cpcm_corr` |
| `7916dc2b20c59b6a21e0925cb81613523643ac6b15631835df97aaac69dab485` | 1795600 | `ino3-r12/workspace/nodes/n0-sp-svp-pbe0/geom-pme3-neutral-min_sp_sp_cpcm_acetonitrile.densities` |
| `c27b436256343a7ccb598c194143d721c4539aa992dd851cb8f8b65c3407de2b` | 1838 | `ino3-r12/workspace/nodes/n0-sp-svp-pbe0/geom-pme3-neutral-min_sp_sp_cpcm_acetonitrile.densitiesinfo` |
| `98366df0fc2dafce8980f8520bbf844ce585c82cce28f92e7ffb0ff5a28d5e80` | 2201668 | `ino3-r12/workspace/nodes/n0-sp-svp-pbe0/geom-pme3-neutral-min_sp_sp_cpcm_acetonitrile.gbw` |
| `e5b142db64f8ef9acd028e3908335bf872d0aa0f4a628c1ae68d7ddbfe2cf96e` | 519727 | `ino3-r12/workspace/nodes/n0-sp-svp-tpssh/geom-pme3-neutral-min_sp_sp_cpcm_acetonitrile.cpcm` |
| `63b08e6550553b0b32f64d17ce8f523b3768885395a2c5837b0fd1ac4878c66a` | 52313 | `ino3-r12/workspace/nodes/n0-sp-svp-tpssh/geom-pme3-neutral-min_sp_sp_cpcm_acetonitrile.cpcm_corr` |
| `e407f5098fe75527b542ca917391218aebb18abbdb49965d3446bc69dd999a07` | 1795600 | `ino3-r12/workspace/nodes/n0-sp-svp-tpssh/geom-pme3-neutral-min_sp_sp_cpcm_acetonitrile.densities` |
| `2e9675862e760d2d41beb28d41416da6f9cd792e4141d7e93669483deed4169a` | 1838 | `ino3-r12/workspace/nodes/n0-sp-svp-tpssh/geom-pme3-neutral-min_sp_sp_cpcm_acetonitrile.densitiesinfo` |
| `7b884a09fb646dc072435617ad326d0b579edf54ecc4604b5d2478512b1afb8d` | 2201668 | `ino3-r12/workspace/nodes/n0-sp-svp-tpssh/geom-pme3-neutral-min_sp_sp_cpcm_acetonitrile.gbw` |
| `335aaacfb73269525aa6e3138dd2e3ade613f573a6917e4dbec2aa5af8b4bb35` | 519727 | `ino3-r12/workspace/nodes/n0-sp-tzvp/geom-pme3-neutral-min_sp_sp_cpcm_acetonitrile.cpcm` |
| `f99e29e6fd07922dd60c6066a2a1ec1ca75883aa7690f0c21aff4e2eb7ceb167` | 52313 | `ino3-r12/workspace/nodes/n0-sp-tzvp/geom-pme3-neutral-min_sp_sp_cpcm_acetonitrile.cpcm_corr` |
| `5a072f26d85654822b9ebcc5ce1e33c6c6819a1491ce36271da452f1fd370909` | 5475600 | `ino3-r12/workspace/nodes/n0-sp-tzvp/geom-pme3-neutral-min_sp_sp_cpcm_acetonitrile.densities` |
| `650d3a067c2894e6adf33e9901d7d78c53770b3691b65a5f2cf9d508b1cf039f` | 1838 | `ino3-r12/workspace/nodes/n0-sp-tzvp/geom-pme3-neutral-min_sp_sp_cpcm_acetonitrile.densitiesinfo` |
| `b7633aa60ef54317829de3e36a74d9ad9a28698a62f711578580f3b3e2134d56` | 4107508 | `ino3-r12/workspace/nodes/n0-sp-tzvp/geom-pme3-neutral-min_sp_sp_cpcm_acetonitrile.gbw` |
| `5a3167b362a8658ee229fe432ddbaa23ca6aeef03b6a53c741aa0b7c649276b4` | 490000 | `po3-r19/workspace/nodes/opt-benzyl-azide/benzyl-azide_opt_opt.densities` |
| `68adab985aa85540f71edad2a2d4614962d7b9c7107574bce9b756924ed0ff50` | 1838 | `po3-r19/workspace/nodes/opt-benzyl-azide/benzyl-azide_opt_opt.densitiesinfo` |
| `9000edf8532ece3c5182d293c5b3522aafa907e5b3ff9243262cb3a5913eb4db` | 1399060 | `po3-r19/workspace/nodes/opt-benzyl-azide/benzyl-azide_opt_opt.gbw` |
| `7d8ec974f448b1c67a56ad56b2a192e84d011bbda0cc87b1310a7597d755562d` | 11776 | `po3-r19/workspace/nodes/opt-benzyl-azide/benzyl-azide_opt_opt.opt` |
| `9c223516f6cae120a7c6cf69d445140ad99c349cba9f09a1aa4e8ce065bebcfb` | 384400 | `po3-r19/workspace/nodes/opt-tfm-ynoate/trifluoromethyl-ynoate_opt_opt.densities` |
| `319f232f95a5d3985d90d9a9a72c004293fb8bc151d8372663e8257e532acd25` | 1838 | `po3-r19/workspace/nodes/opt-tfm-ynoate/trifluoromethyl-ynoate_opt_opt.densitiesinfo` |
| `92476da91284d9b2a04e835a315603043bfc9ea8b924a1414c2096008959d5a7` | 1325172 | `po3-r19/workspace/nodes/opt-tfm-ynoate/trifluoromethyl-ynoate_opt_opt.gbw` |
| `c3641bed7b819497fb8a98f5b417717c5735adc17add6c13fbe275a2b044a556` | 7556 | `po3-r19/workspace/nodes/opt-tfm-ynoate/trifluoromethyl-ynoate_opt_opt.opt` |
| `0ab6e603f7e77d1b1c47b4e70ef52cde22d5b7403f0912dc7cffb3f11fbdbf26` | 2177188 | `po3-r19/workspace/nodes/scan-esterc4-path/ts-guess-esterc4-complex_scan_scan.001.gbw` |
| `8efffaf29d92606ef219609ffc290b867c7dbcf44b140addaa0e4be7ad029a8a` | 2177188 | `po3-r19/workspace/nodes/scan-esterc4-path/ts-guess-esterc4-complex_scan_scan.002.gbw` |
| `7bc5ffb5e3ef5be05e00b093be8717ea6c2c65cf50713916d49b2c511c7491a7` | 2177188 | `po3-r19/workspace/nodes/scan-esterc4-path/ts-guess-esterc4-complex_scan_scan.003.gbw` |
| `f66a6404bd75b25cc9dd22d24a2cb36a0733af26e96f773358a8094c09f07b97` | 2177188 | `po3-r19/workspace/nodes/scan-esterc4-path/ts-guess-esterc4-complex_scan_scan.004.gbw` |
| `452cdae0f6eda596da608474d8446ea8588a50a08dfad64f43e0a93c47933461` | 2177188 | `po3-r19/workspace/nodes/scan-esterc4-path/ts-guess-esterc4-complex_scan_scan.005.gbw` |
| `77424edbeaf40755f3e698e6d30638a21f0ddd7ccd396b0f69c37f5d33ebfa1e` | 2177188 | `po3-r19/workspace/nodes/scan-esterc4-path/ts-guess-esterc4-complex_scan_scan.006.gbw` |
| `4f5d1f2f1b620a0e5219237889f2347023a6148738e2871387abd2b7ceb331c5` | 2177188 | `po3-r19/workspace/nodes/scan-esterc4-path/ts-guess-esterc4-complex_scan_scan.007.gbw` |
| `f894388745750f4e439a6b740e5ac3392dc1c5ca663deb6b2e7e9d0d53af944c` | 1742400 | `po3-r19/workspace/nodes/scan-esterc4-path/ts-guess-esterc4-complex_scan_scan.densities` |
| `c5f63b91e9f7b4ad794b286f9960d8504c8837d5e0941196a8e2d6fc6c15b9d1` | 1838 | `po3-r19/workspace/nodes/scan-esterc4-path/ts-guess-esterc4-complex_scan_scan.densitiesinfo` |
| `4f5d1f2f1b620a0e5219237889f2347023a6148738e2871387abd2b7ceb331c5` | 2177188 | `po3-r19/workspace/nodes/scan-esterc4-path/ts-guess-esterc4-complex_scan_scan.gbw` |
| `9629447eb00257c0c85b3e1d1d134006f612279dc4c90b9e88d7dd0e8df9797b` | 33596 | `po3-r19/workspace/nodes/scan-esterc4-path/ts-guess-esterc4-complex_scan_scan.opt` |
| `f2cc9ec983c7c7c0fd50baf489e78696d035ed727cf32ddd7b0b6cf0ea143a05` | 2177188 | `po3-r19/workspace/nodes/scan-esterc5-path/ts-guess-esterc4-complex_scan_scan.001.gbw` |
| `5039a83086c1fac831f0023f7bf54351b5596677c75544707c5a46e4ac9c6ce3` | 2177188 | `po3-r19/workspace/nodes/scan-esterc5-path/ts-guess-esterc4-complex_scan_scan.002.gbw` |
| `a61d341e3a35127e8ad2bb3b123f26e7d70259dc34c2ad70396b61cb4dca3a4b` | 2177188 | `po3-r19/workspace/nodes/scan-esterc5-path/ts-guess-esterc4-complex_scan_scan.003.gbw` |
| `350ade13ac2eb375d4ee73545b3648b7498a2952f09dd9180ea59b4058b619dd` | 2177188 | `po3-r19/workspace/nodes/scan-esterc5-path/ts-guess-esterc4-complex_scan_scan.004.gbw` |
| `85f2c392dc4015f69299a3c4bf665aac05f454ca3f725101bbbcdac46d3dda92` | 2177188 | `po3-r19/workspace/nodes/scan-esterc5-path/ts-guess-esterc4-complex_scan_scan.005.gbw` |
| `57deee4cde3c1c9486cb09e92c666edaac8c23bf58bf49ec52b303ce3cdbad38` | 2177188 | `po3-r19/workspace/nodes/scan-esterc5-path/ts-guess-esterc4-complex_scan_scan.006.gbw` |
| `e6502f8b2d5aaaf7385a36e2fa9f210abacb49bdeb273331cd8a6542a70d21e6` | 2177188 | `po3-r19/workspace/nodes/scan-esterc5-path/ts-guess-esterc4-complex_scan_scan.007.gbw` |
| `83c6e72d97e857842e881e778b11bdb35c8e86cf892018a458894a9a7bf047d9` | 1742400 | `po3-r19/workspace/nodes/scan-esterc5-path/ts-guess-esterc4-complex_scan_scan.densities` |
| `2d4130f140c258876d6dc685a8d1066fe73480eb94579ddee973b603a683f569` | 1838 | `po3-r19/workspace/nodes/scan-esterc5-path/ts-guess-esterc4-complex_scan_scan.densitiesinfo` |
| `e6502f8b2d5aaaf7385a36e2fa9f210abacb49bdeb273331cd8a6542a70d21e6` | 2177188 | `po3-r19/workspace/nodes/scan-esterc5-path/ts-guess-esterc4-complex_scan_scan.gbw` |
| `3811cdd653b55b73aad04cbadecb37d3540d75e2c56687b75cc7b24c549692f4` | 31036 | `po3-r19/workspace/nodes/scan-esterc5-path/ts-guess-esterc4-complex_scan_scan.opt` |
| `2b4a10a72ad34fa980bd908f31229d141c755029e7390635a2b6bc2dff240ba0` | 7398400 | `po3-r19/workspace/nodes/sp-c4-b3lyp-tzvp/geom-c4-saddle-sp_sp_sp_gas_phase.densities` |
| `4207d8be05bb81d8b06431511e85174a39dd1e9b7674fe5c76d3f58dcb0047cf` | 1838 | `po3-r19/workspace/nodes/sp-c4-b3lyp-tzvp/geom-c4-saddle-sp_sp_sp_gas_phase.densitiesinfo` |
| `45bd7b8ad9de29f7d1bb3b2976885cf8b92246ccb7dc6a41a4b8afb3686b2730` | 5088356 | `po3-r19/workspace/nodes/sp-c4-b3lyp-tzvp/geom-c4-saddle-sp_sp_sp_gas_phase.gbw` |
| `f12246e0accbd662b5d8753bed1f4937e612f34bc194cf75e96679e6e11fec98` | 7398400 | `po3-r19/workspace/nodes/sp-c4-dsd-tzvp/geom-c4-saddle-sp_sp_sp_gas_phase.densities` |
| `a7ade7dd0560595749c96cbd9edc5382c13763b272a804f016bbdbbd5ad604e2` | 1838 | `po3-r19/workspace/nodes/sp-c4-dsd-tzvp/geom-c4-saddle-sp_sp_sp_gas_phase.densitiesinfo` |
| `c34dd9c66b1b8849d0c9a1a9f96a984c78bbb162ef71e69dcded4731b4009604` | 7721716 | `po3-r19/workspace/nodes/sp-c4-dsd-tzvp/geom-c4-saddle-sp_sp_sp_gas_phase.gbw` |
| `621d17ffc76e5ef438f3cf8230ac8af449d9ed82955b44ba8aa1414392c78670` | 7398400 | `po3-r19/workspace/nodes/sp-c5-b3lyp-tzvp/geom-c5-saddle-sp_sp_sp_gas_phase.densities` |
| `a3fb3bac22d40e1ff9d1d922fb27fd008bbdcbb765f31c99f68f5c3252c2d379` | 1838 | `po3-r19/workspace/nodes/sp-c5-b3lyp-tzvp/geom-c5-saddle-sp_sp_sp_gas_phase.densitiesinfo` |
| `f6d51664916b2cf2ff33cacdfb8facdf0379256d29ee98c8292a47ff6397dca4` | 5088356 | `po3-r19/workspace/nodes/sp-c5-b3lyp-tzvp/geom-c5-saddle-sp_sp_sp_gas_phase.gbw` |
| `fe9e7fafbd289fced3c866641d6c5749fcfb63fe27b296a226276dbc407c7b5f` | 7398400 | `po3-r19/workspace/nodes/sp-c5-dsd-tzvp/geom-c5-saddle-sp_sp_sp_gas_phase.densities` |
| `3626273446d22d4a879386c5849b159ebdbbbea60e2b16380d4fc08807775b5c` | 1838 | `po3-r19/workspace/nodes/sp-c5-dsd-tzvp/geom-c5-saddle-sp_sp_sp_gas_phase.densitiesinfo` |
| `db0bfd19fff6443f5a802b33dba5fda43e5b05e4f1a394b51a606b515a0b3067` | 7721716 | `po3-r19/workspace/nodes/sp-c5-dsd-tzvp/geom-c5-saddle-sp_sp_sp_gas_phase.gbw` |
| `630b6767fed1a05a74188bc45da5dc941a28f0424c5e67cf1a0311819a837560` | 1742400 | `po3-r19/workspace/nodes/ts-esterc4/presaddle-esterc4_optts_optts.densities` |
| `456e42637f0da69dc7ec343a1ea08f2ce541de6aec90710930d89660e0043ab2` | 1838 | `po3-r19/workspace/nodes/ts-esterc4/presaddle-esterc4_optts_optts.densitiesinfo` |
| `1af2a7e68a22fade95a927e4800fec315e0fbdc6ec685f14839e80a41429c90f` | 2177156 | `po3-r19/workspace/nodes/ts-esterc4/presaddle-esterc4_optts_optts.gbw` |
| `05b815176da29370fe15ea546bce65676d73234b13f1c81d3985ef4d1df8b01d` | 22628 | `po3-r19/workspace/nodes/ts-esterc4/presaddle-esterc4_optts_optts.opt` |
| `bc4ccaded4438b358ce06853107d6e31ce2b863a84c0e40a6cd5f3b5de000c3e` | 1742400 | `po3-r19/workspace/nodes/ts-esterc4-restart/geom-c4-saddle-reached_optts_optts.densities` |
| `13deb617bc6421a83b6a549782b40f1e3d1886a7cb257669e5c7e32a4478ad12` | 1838 | `po3-r19/workspace/nodes/ts-esterc4-restart/geom-c4-saddle-reached_optts_optts.densitiesinfo` |
| `6c8f334a17c56c2302ee6a09eda778d138c2770a8959bdcfbfb3db1f724ff73f` | 2177156 | `po3-r19/workspace/nodes/ts-esterc4-restart/geom-c4-saddle-reached_optts_optts.gbw` |
| `3884da96bbb5fdca056e24f0eb21e1da90a702e02ee279f5742496a11bbd8fdb` | 9236 | `po3-r19/workspace/nodes/ts-esterc4-restart/geom-c4-saddle-reached_optts_optts.opt` |
| `5369fa97e4406b449b6f1d4d1cf474e5d5d861f657c56dd057635ab737b91aaa` | 1742400 | `po3-r19/workspace/nodes/ts-esterc5/presaddle-esterc5_optts_optts.densities` |
| `334a0631e17acbff11afa88fd17dc6b48492f3a6f8bc3f9ab981f692cdc7cd25` | 1838 | `po3-r19/workspace/nodes/ts-esterc5/presaddle-esterc5_optts_optts.densitiesinfo` |
| `fdf11479ba291420a1e429141fcc51cf5f59eb645e9483e9bcb4fca270827da7` | 2177156 | `po3-r19/workspace/nodes/ts-esterc5/presaddle-esterc5_optts_optts.gbw` |
| `f4e1be506c599804651fbda1ec1e0b075cae3d761449de7005a51ecac668fb74` | 23484 | `po3-r19/workspace/nodes/ts-esterc5/presaddle-esterc5_optts_optts.opt` |
