data:
	sh run_phonopy_projector_for_MAPIBrCl.sh

plots: mode_decompositions.gpt
	gnuplot mode_decompositions.gpt

montage: MAPbI_mode.png  MAPbBr_mode.png  MAPbCl_mode.png
	montage MAPbI_mode.png  MAPbBr_mode.png  MAPbCl_mode.png  -geometry 1024x1024 MAPbI_Br_Cl_mode_montage_horizontal.png
	montage -tile 1x MAPbI_mode.png  MAPbBr_mode.png  MAPbCl_mode.png  -geometry 1024x1024 MAPbI_Br_Cl_mode_montage_vertical.png




tightmontage: MAPbI_mode.png  MAPbBr_mode.png  MAPbCl_mode.png
	montage MAPbI_mode.png -geometry 1024x1024+0x0 MAPbBr_mode.png  -crop 900x1024+100x0 MAPbCl_mode.png -crop 900x1024+100x0  -geometry 900x1024 MAPbI_Br_Cl_mode_tightmontage.png

clean:
	rm */anim_??.xyz

