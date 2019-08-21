# (1) register log(a^x) = x*log(a)
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

# (2) register exp(x)*exp(y) -> exp(x+y)
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

# (3) register (a^{x})^{b} = (a^{b})^{x}
src_nds = Dict{Int, Template_Node}(1 => Template_Node(:op, :^),
                                   2 => Template_Node(:op, :^),
                                   3 => Template_Node(:num, :b),
                                   4 => Template_Node(:num, :a),
                                   5 => Template_Node(:expr, :x))
src_dag = [5 => 2, 4 => 2, 3 => 1, 2 => 1]
src = Template_Graph(src_nds, src_dag)
dest_nds = Dict{Int, Template_Node}(1 => Template_Node(:op, :^),
                                    2 => Template_Node(:op, :^),
                                    3 => Template_Node(:expr, :x),
                                    4 => Template_Node(:num, :a),
                                    5 => Template_Node(:num, :b))
dest_dag = [5 => 2, 4 => 2, 3 => 1, 2 => 1]
dest = Template_Graph(dest_nds, dest_dag)
register_substitution!(src, dest)

# (4) register (x^{a})^{b} = x^{(ab)}
src_nds = Dict{Int, Template_Node}(1 => Template_Node(:op, :^),
                                   2 => Template_Node(:op, :^),
                                   3 => Template_Node(:num, :a),
                                   4 => Template_Node(:expr, :x),
                                   5 => Template_Node(:num, :b))
src_dag = [5 => 2, 4 => 2, 3 => 1, 2 => 1]
src = Template_Graph(src_nds, src_dag)
dest_nds = Dict{Int, Template_Node}(1 => Template_Node(:op, :^),
                                    2 => Template_Node(:expr, :x),
                                    3 => Template_Node(:op, :*),
                                    4 => Template_Node(:num, :a),
                                    5 => Template_Node(:num, :b))
dest_dag = [5 => 3, 4 => 3, 3 => 1, 2 => 1]
dest = Template_Graph(dest_nds, dest_dag)
register_substitution!(src, dest)

# (5) register a^{\log(x)} = x^{\log(a)}
src_nds = Dict{Int, Template_Node}(1 => Template_Node(:op, :^),
                                   2 => Template_Node(:num, :a),
                                   3 => Template_Node(:op, :log),
                                   4 => Template_Node(:expr, :x))
src_dag = [4 => 3, 3 => 1, 2 => 1]
src = Template_Graph(src_nds, src_dag)
dest_nds = Dict{Int, Template_Node}(1 => Template_Node(:op, :^),
                                    2 => Template_Node(:expr, :x),
                                    3 => Template_Node(:op, :log),
                                    4 => Template_Node(:num, :a))
dest_dag = [4 => 3, 3 => 1, 2 => 1]
dest = Template_Graph(dest_nds, dest_dag)
register_substitution!(src, dest)

# (6) register \log(xy) = \log(x) + \log(y) DONEish
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

# (7) register \log(x/y) = \log(x) - \log(y)
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
