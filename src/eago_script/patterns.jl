# 1 register log(a^x) = x*log(a) DONEish
println("pattern #1")
src_nds = Dict{Int, Template_Node}(1 => Template_Node(:op, :log),
                                   2 => Template_Node(:op, :^),
                                   3 => Template_Node(:num, :a; check = a -> a >= 0.0),
                                   4 => Template_Node(:expr, :x))
src_dag = [4 => 2, 3 => 2, 2 => 1]
src = Template_Graph(src_nds, src_dag)
dest_nds = Dict{Int, Template_Node}(1 => Template_Node(:op, :*),
                                    2 => Template_Node(:expr, :x),
                                    3 => Template_Node(:op, :log),
                                    4 => Template_Node(:num, :a))
dest_dag = [4 => 3, 3 => 1, 2 => 1]
dest = Template_Graph(dest_nds, dest_dag)
register_substitution!(src, dest)


# 2 register exp(x)*exp(y) DONEish
println("pattern #2")
src_nds = Dict{Int, Template_Node}(1 => Template_Node(:op, :*),
                                   2 => Template_Node(:op, :exp),
                                   3 => Template_Node(:op, :exp),
                                   4 => Template_Node(:expr, :x),
                                   5 => Template_Node(:expr, :y))
src_dag = [5 => 3, 4 => 2, 3 => 1, 2 => 1]
src = Template_Graph(src_nds, src_dag)
dest_nds = Dict{Int, Template_Node}(1 => Template_Node(:op, :exp),
                                    2 => Template_Node(:op, :+),
                                    3 => Template_Node(:expr, :x),
                                    4 => Template_Node(:expr, :y))
dest_dag = [4 => 2, 3 => 2, 2 => 1]
dest = Template_Graph(dest_nds, dest_dag)
register_substitution!(src, dest)


# 3 (a^{x})^{b} = (a^{b})^{x} DONEish
println("pattern #3")
src_nds = Dict{Int, Template_Node}(1 => Template_Node(:op, :^),
                                   2 => Template_Node(:op, :^),
                                   3 => Template_Node(:num, :b),
                                   4 => Template_Node(:num, :a),
                                   5 => Template_Node(:expr, :x))
src_dag = [5 => 2, 4 => 2, 2 => 1, 3 => 1]
dest_nds = Dict{Int, Template_Node}(1 => Template_Node(:op, :^),
                                    2 => Template_Node(:op, :^),
                                    3 => Template_Node(:expr, :x),
                                    5 => Template_Node(:num, :b),
                                    4 => Template_Node(:num, :a))
dest_dag = [5 => 2, 4 => 2, 3 => 1, 2 => 1]
dest = Template_Graph(dest_nds, dest_dag)
register_substitution!(src, dest)


# 4 (x^{a})^{b} = x^{(ab)} DONEish
println("pattern #4")
src_nds = Dict{Int, Template_Node}(1 => Template_Node(:op, :^),
                                   2 => Template_Node(:num, :b),
                                   3 => Template_Node(:op, :^),
                                   5 => Template_Node(:num, :a),
                                   4 => Template_Node(:expr, :x))
src_dag = [5 => 3, 4 => 3, 3 => 1, 2 => 1]
src = Template_Graph(src_nds, src_dag)
dest_nds = Dict{Int, Template_Node}(1 => Template_Node(:op, :^),
                                    3 => Template_Node(:op, :*),
                                    4 => Template_Node(:num, :a),
                                    5 => Template_Node(:num, :b),
                                    2 => Template_Node(:expr, :x))
dest_dag = [5 => 3, 4 => 3, 3 => 1, 2 => 1]
dest = Template_Graph(dest_nds, dest_dag)
register_substitution!(src, dest)


# 5 a^{\log(x)} = x^{\log(a)} DONEish
println("pattern #5")
src_nds = Dict{Int, Template_Node}(1 => Template_Node(:op, :^),
                                   2 => Template_Node(:op, :log),
                                   3 => Template_Node(:expr, :x),
                                   4 => Template_Node(:num, :a))
src_dag = [4 => 1, 3 => 2, 2 => 1]
src = Template_Graph(src_nds, src_dag)
dest_nds = Dict{Int, Template_Node}(1 => Template_Node(:op, :^),
                                    2 => Template_Node(:expr, :x),
                                    3 => Template_Node(:op, :log),
                                    4 => Template_Node(:num, :a))
dest_dag = [4 => 3, 3 => 1, 2 => 1]
dest = Template_Graph(dest_nds, dest_dag)
register_substitution!(src, dest)


# 6 \log(xy) = \log(x) + \log(y) DONEish
println("pattern #6")
src_nds = Dict{Int, Template_Node}(1 => Template_Node(:op, :log),
                                   2 => Template_Node(:op, :*),
                                   3 => Template_Node(:expr, :x),
                                   4 => Template_Node(:expr, :y))
src_dag = [4 => 2, 3 => 2, 2 => 1]
src = Template_Graph(src_nds, src_dag)
dest_nds = Dict{Int, Template_Node}(1 => Template_Node(:op, :+),
                                    2 => Template_Node(:op, :log),
                                    3 => Template_Node(:op, :log),
                                    4 => Template_Node(:expr, :x),
                                    5 => Template_Node(:expr, :y))
dest_dag = [5 => 3, 4 => 2, 3 => 1, 2 => 1]
dest = Template_Graph(dest_nds, dest_dag)
register_substitution!(src, dest)

#=
# \log(x/y) = \log(x) - \log(y)  DONEish
println("pattern #7")
src_nds = Dict{Int, Template_Node}(1 => Template_Node(:op, :log),
                                   2 => Template_Node(:op, :/),
                                   3 => Template_Node(:expr, :x),
                                   4 => Template_Node(:expr, :y))
src_dag = [4 => 2, 3 => 2, 2 => 1]
src = Template_Graph(src_nds, src_dag)
dest_nds = Dict{Int, Template_Node}(1 => Template_Node(:op, :-),
                                    2 => Template_Node(:op, :log),
                                    3 => Template_Node(:op, :log),
                                    4 => Template_Node(:expr, :x),
                                    5 => Template_Node(:expr, :y))
dest_dag = [5 => 3, 4 => 2, 3 => 1, 2 => 1]
dest = Template_Graph(dest_nds, dest_dag)
register_substitution!(src, dest)
=#
