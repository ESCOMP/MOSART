module RtmMctMod

 use m_MCTWorld         , only: mct_world_init    => init

 use m_AttrVect         , only: mct_aVect         => AttrVect
 use m_AttrVect         , only: mct_aVect_init    => init
 use m_AttrVect         , only: mct_aVect_clean   => clean
 use m_AttrVect         , only: mct_aVect_zero    => zero
 use m_AttrVect         , only: mct_aVect_lsize   => lsize
 use m_AttrVect         , only: mct_aVect_indexIA => indexIA
 use m_AttrVect         , only: mct_aVect_indexRA => indexRA

 use m_AttrVectComms    , only: mct_aVect_scatter => scatter
 use m_AttrVectComms    , only: mct_aVect_gather  => gather

 use m_MatAttrVectMul   , only: mct_sMat_avMult   => sMatAvMult

 use m_GlobalSegMap     , only: mct_gsMap         => GlobalSegMap
 use m_GlobalSegMap     , only: mct_gsMap_init    => init
 use m_GlobalSegMap     , only: mct_gsMap_clean   => clean

 use m_SparseMatrix     , only: mct_sMat          => SparseMatrix
 use m_SparseMatrix     , only: mct_sMat_Init     => init
 use m_SparseMatrix     , only: mct_sMat_Clean    => clean
 use m_SparseMatrix     , only: mct_sMat_indexIA  => indexIA
 use m_SparseMatrix     , only: mct_sMat_indexRA  => indexRA
 use m_SparseMatrix     , only: mct_sMat_GNumEl   => GlobalNumElements

 use m_SparseMatrixPlus , only: mct_sMatP         => SparseMatrixPlus
 use m_SparseMatrixPlus , only: mct_sMatP_Init    => init
 use m_SparseMatrixPlus , only: mct_sMatP_clean   => clean

end module RtmMctMod
