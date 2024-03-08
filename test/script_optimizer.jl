function check_node(nd::EAGO.Script.NodeInfo, type::MOINL.NodeType , indx::Int, child::Vector{Int})
    flag = true
    if nd.nodetype != type
        flag = false
        return flag
    end
    if nd.index != indx
        flag = false
        return flag
    end
    for i in 1:length(child)
        if nd.children[i] != child[i]
            flag = false
            return flag
        end
    end
    return flag
end

@testset "Test Trace Script" begin
    op_plus = EAGO.Script.multivariate_operator_to_id[:+]
    op_mul = EAGO.Script.multivariate_operator_to_id[:*]
    op_sin = EAGO.Script.univariate_operator_to_id[:sin]
    op_abs = EAGO.Script.univariate_operator_to_id[:abs]
    function f1(x)
        return sin(3.0*x[1]) + x[2]
    end
    tape1 = EAGO.Script.trace_script(f1,2)
    @test check_node(tape1.nd[1], MOINL.NODE_VARIABLE, 1, [-1])
    @test check_node(tape1.nd[2], MOINL.NODE_VARIABLE, 2, [-1])
    @test check_node(tape1.nd[3], MOINL.NODE_VALUE, 1, [-2])
    @test check_node(tape1.nd[4], MOINL.NODE_CALL_MULTIVARIATE, op_mul, [3,1])
    @test check_node(tape1.nd[5], MOINL.NODE_CALL_UNIVARIATE, op_sin, [4])
    @test check_node(tape1.nd[6], MOINL.NODE_CALL_MULTIVARIATE, op_plus, [5,2])
    @test tape1.const_values[1] == 3.0
    @test tape1.num_valued[1] == false
    @test tape1.num_valued[2] == false
    @test tape1.num_valued[3] == true
    @test tape1.num_valued[4] == false
    @test tape1.num_valued[5] == false
    @test tape1.num_valued[6] == false
    @test tape1.set_trace_count == 6
    @test tape1.const_count == 1

    function f2(x)
        z = abs(x[1])
        y = 1.3
        for i in 1:2
            y += i*z
        end
        return sin(3.0*x[2]) + y
    end
    tape2 = EAGO.Script.trace_script(f2,2)
    @test check_node(tape2.nd[1], MOINL.NODE_VARIABLE, 1, [-1])
    @test check_node(tape2.nd[2], MOINL.NODE_VARIABLE, 2, [-1])
    @test check_node(tape2.nd[3], MOINL.NODE_CALL_UNIVARIATE, op_abs, [1])
    @test check_node(tape2.nd[4], MOINL.NODE_VALUE, 1, [-2])
    @test check_node(tape2.nd[5], MOINL.NODE_CALL_MULTIVARIATE, op_mul, [4,3])
    @test check_node(tape2.nd[6], MOINL.NODE_VALUE, 2, [-2])
    @test check_node(tape2.nd[7], MOINL.NODE_CALL_MULTIVARIATE, op_plus, [6,5])
    @test check_node(tape2.nd[8], MOINL.NODE_VALUE, 3, [-2])
    @test check_node(tape2.nd[9], MOINL.NODE_CALL_MULTIVARIATE, op_mul, [8,3])
    @test check_node(tape2.nd[10], MOINL.NODE_CALL_MULTIVARIATE, op_plus, [7,9])
    @test check_node(tape2.nd[11], MOINL.NODE_VALUE, 4, [-2])
    @test check_node(tape2.nd[12], MOINL.NODE_CALL_MULTIVARIATE, op_mul, [11,2])
    @test check_node(tape2.nd[13], MOINL.NODE_CALL_UNIVARIATE, op_sin, [12])
    @test check_node(tape2.nd[14], MOINL.NODE_CALL_MULTIVARIATE, op_plus, [13,10])
    @test tape2.const_values[1] == 1.0
    @test tape2.const_values[2] == 1.3
    @test tape2.const_values[3] == 2.0
    @test tape2.const_values[4] == 3.0
    @test tape2.num_valued[1] == false
    @test tape2.num_valued[2] == false
    @test tape2.num_valued[3] == false
    @test tape2.num_valued[4] == true
    @test tape2.num_valued[5] == false
    @test tape2.num_valued[6] == true
    @test tape2.num_valued[7] == false
    @test tape2.num_valued[8] == true
    @test tape2.num_valued[9] == false
    @test tape2.num_valued[10] == false
    @test tape2.num_valued[11] == true
    @test tape2.num_valued[12] == false
    @test tape2.num_valued[13] == false
    @test tape2.num_valued[14] == false
    @test tape2.set_trace_count == 14
    @test tape2.const_count == 4

    #=
    nds = EAGO.Script.tape_to_list(tape2)
    @test nds[1].index == 1
    @test nds[1].children[1] == -1
    @test nds[13].index == 15
    @test nds[13].children[1] == 12
    =#

    function f3(x)
        z = abs(x[1])::Float64
        return z
    end
    tape3 = EAGO.Script.trace_script(f3,2)
    #=
    @test tape3.const_count == 0
    @test tape3.set_trace_count == 3
    @test check_node(tape3.nd[1], MOINL.NODE_VARIABLE, 1, [-1])
    @test check_node(tape3.nd[2], MOINL.NODE_VARIABLE, 2, [-1])
    @test check_node(tape3.nd[3], MOINL.NODE_CALL_UNIVARIATE, 3, [1])
    =#
end
