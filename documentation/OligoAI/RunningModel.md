To run the Oligo AI, we installed the envrionment.yml that was given.

Since the GPU available was RTX 5090, we uninstalled the existing flash-attention library and used a new beta version instead 
https://github.com/Dao-AILab/flash-attention/releases/tag/fa4-v4.0.0.beta4
that supports the Blackwell architecture. 

To fix minor issues caused by this move we forked the OligoAI repository and made backwards compatible changes
https://github.com/RedPenguin100/OligoAI

