using SnoopCompileCore, SnoopCompile
trees = invalidation_trees(@snoopr using EAGO)
methinvs = trees[end]
root = methinvs.backedges[end]
ascend(root)