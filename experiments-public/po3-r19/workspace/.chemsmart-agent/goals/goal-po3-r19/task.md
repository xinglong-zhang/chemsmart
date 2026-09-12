I am about to run the uncatalysed thermal cycloaddition of benzyl azide onto methyl 4,4,4-trifluorobut-2-ynoate, neat at 80 C. The alkyne takes me four steps and I have 300 mg of it, so I would rather know in advance whether this gives me one regioisomer or a mixture. Both ends of the triple bond are electron-poor and I cannot decide which one controls the orientation: the ester is the conventional activating group and I would expect it to end up at C4, but the trifluoromethyl is polarising the alkyne as well and I do not know whether it reinforces or opposes that. I need the isomer bearing the ester at C4 for the next step. So: which regiochemistry does the thermal reaction favour, and by how much in activation free energy at 353 K? If the two orientations are within a few tenths of a kcal/mol I will not bother with the thermal route and will go straight to a ruthenium catalyst instead.

What would satisfy me: the major regioisomer named unambiguously, plus ddG(activation) between the two orientations in kcal/mol at 353 K. Anything below roughly 0.5 kcal/mol I will read as "no thermal condition will fix this".

Structures in the workspace (hand-built starts from SMILES, RDKit plus MMFF, nothing optimised; all neutral singlets):
  benzyl-azide.xyz             [N-]=[N+]=NCc1ccccc1
  trifluoromethyl-ynoate.xyz   COC(=O)C#CC(F)(F)F
  triazole-ester-at-c4.xyz     COC(=O)c1c(C(F)(F)F)nnn1Cc2ccccc2
  triazole-ester-at-c5.xyz     COC(=O)c1nnn(Cc2ccccc2)c1C(F)(F)F