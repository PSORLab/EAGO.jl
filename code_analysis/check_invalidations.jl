using SnoopCompile
invalidations = @snoopr using EAGO

@show length(invalidations)
trees = invalidation_trees(invalidations)
#methinvs = trees[end]
#root = methinvs.backedges[end]
#ascend(root)