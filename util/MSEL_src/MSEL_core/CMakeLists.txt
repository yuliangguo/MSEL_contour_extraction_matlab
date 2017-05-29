
	INCLUDE_DIRECTORIES( ${VXL_GEL_INCLUDE_DIR} )
	INCLUDE_DIRECTORIES( ${VXL_MUL_INCLUDE_DIR} )	
	#MESSAGE("VXL_INCLUDE_DIR is ${VXL_CORE_INCLUDE_DIR}")
	#MESSAGE("VXL_GEL_INCLUDE_DIR is ${VXL_GEL_INCLUDE_DIR}")
	#MESSAGE("VXL_MUL_INCLUDE_DIR is ${VXL_MUL_INCLUDE_DIR}")
	set( MSEL_core_sources
		# From bgld and bnld in vxd
		bgld_eulerspiral_base.h
		bgld_eulerspiral.h			bgld_eulerspiral.cxx
		bgld_arc.h 				bgld_arc.cxx
		bgld_arc_base.h 			bgld_arc_base.cxx
		bgld_arc_algo.h 			bgld_arc_algo.cxx
		bgld_biarc.h 				bgld_biarc.cxx
		bgld_curve_smoothing.h			bgld_curve_smoothing.cxx
		bgld_param_curve.h
		bnld_fresnel.h				bnld_fresnel.cxx

		# From breye1d in vxd
		bpro1_parameters.h			bpro1_parameters.cxx			bpro1_parameters_sptr.h
		bpro1_parameters.hxx		
		bpro1_process.h 			bpro1_process.cxx			bpro1_process_sptr.h
		bpro1_storage.h			 	bpro1_storage.cxx			bpro1_storage_sptr.h
		vidpro1_vsol2D_storage.h		vidpro1_vsol2D_storage.cxx		vidpro1_vsol2D_storage_sptr.h

		# Templates from breye1d in vxd
		bpro1_parameters+bool-.cxx
		bpro1_parameters+bpro_filepath-.cxx
		bpro1_parameters+double-.cxx
		bpro1_parameters+unsigned-.cxx
		vbl_smart_ptr+bpro1_storage-.cxx
		vbl_smart_ptr+bpro1_process-.cxx
		vbl_smart_ptr+bpro1_parameters-.cxx				
		vsl_map_io+vcl_string.vcl_vector+vsol_spatial_object_2d_sptr--.cxx
		vsl_map_io+vcl_string.vcl_vector+double--.cxx
		vsl_binary_loader+bpro1_storage-.cxx



		# dbdet codes from VXL
		dbdet_appearance.h
		dbdet_cem_file_io.h			dbdet_cem_file_io.cxx
		dbdet_CFTG.h
		dbdet_curve_fragment_graph.h
		dbdet_curvelet.h			dbdet_curvelet.cxx
		dbdet_curvelet_map.h			dbdet_curvelet_map.cxx
		dbdet_curve_model.h			dbdet_curve_model.cxx
		dbdet_edgel.h				dbdet_edgel.cxx
		dbdet_edgel_link_graph.h
		dbdet_edgemap.h				dbdet_edgemap.cxx			dbdet_edgemap_sptr.h
		dbdet_edgemap_storage.h			dbdet_edgemap_storage.cxx		dbdet_edgemap_storage_sptr.h
		dbdet_EHT.h
		dbdet_EMD.h				dbdet_EMD.cxx
		dbdet_hyp_tree_graph.h
		dbdet_load_edg.h			dbdet_load_edg.cxx
		dbdet_load_edg_process.h		dbdet_load_edg_process.cxx
		dbdet_postprocess_contours.h		dbdet_postprocess_contours.cxx
		dbdet_save_cem_process.h		dbdet_save_cem_process.cxx
		dbdet_sel_base_CFTG_algo.cxx
		dbdet_sel_base.h			dbdet_sel_base.cxx

		dbdet_sel.cxx
		dbdet_sel+dbdet_CC_curve_model_3d-.cxx
		dbdet_sel+dbdet_CC_curve_model-.cxx
		dbdet_sel+dbdet_CC_curve_model_new-.cxx
		dbdet_sel+dbdet_CC_curve_model_perturbed-.cxx
		dbdet_sel+dbdet_ES_curve_model-.cxx
		dbdet_sel+dbdet_ES_curve_model_perturbed-.cxx
		dbdet_sel+dbdet_linear_curve_model-.cxx
		dbdet_sel+dbdet_simple_linear_curve_model-.cxx
		dbdet_sel.h									dbdet_sel_sptr.h
		dbdet_sel_process.h			dbdet_sel_process.cxx
		dbdet_sel_storage.h			dbdet_sel_storage.cxx			dbdet_sel_storage_sptr.h
		dbdet_sel_utils.h

		# Templates from dbdet in vxd
		vbl_smart_ptr+dbdet_edgemap-.cxx
		vbl_smart_ptr+dbdet_sel_base-.cxx

	)
	add_library (MSEL_core ${MSEL_core_sources} )
	target_link_libraries (MSEL_core pdf1d vdgl) 


