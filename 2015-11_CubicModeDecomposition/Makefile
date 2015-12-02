data:
	sh run_phonopy_projector_for_MAPIBrCl.sh

plots: mode_decompositions.gpt
	gnuplot mode_decompositions.gpt

montage: MAPbI_mode.png  MAPbBr_mode.png  MAPbCl_mode.png
	montage MAPbI_mode.png  MAPbBr_mode.png  MAPbCl_mode.png  -geometry 1024x1024 MAPbI_Br_Cl_mode_montage.png

clean:
	rm */anim_??.xyz

