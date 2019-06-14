@testset "Test Base Implicit Routines (Out-of-place)" begin
    EAGO.set_mc_differentiability!(0)
    opts1 = mc_opts(0.5,1,:Dense,:Newton,1,1,1E-10)

    f(x,p) = x[1]*p[1]+p[1]
    g(x,p) = [x[1]*p[1]+p[1];
              x[1]*p[1]+2*p[1]]
    function h1(x,p)
        t1 = x[1]^2
        t2 = x[1]*p[1]
        t3 = 4.0
        t4 = t1 + t2
        t5 = t4 + t3
        return [t5]
    end

    hj1(x,p) = [2*x[1]+p[1]]

    P = [Interval(6.0,9.0)]
    X = [Interval(-0.78,-0.4)]
    p = [7.5]
    pmid = mid.(P)

    xIntv1 = Interval(1.0,3.0)
    xIBox = SVector{1,Interval{Float64}}([xIntv1])
    mBox = mid.(xIBox)
    np = 1
    szero = @SVector zeros(np)
    sone = @SVector ones(np)
    p_mc = [MC{np}(p[i],p[i],Interval(P[i].lo,P[i].hi),sone,sone,false) for i=1:np]

    param = gen_expansion_params(h1,hj1,X,P,pmid,opts1)

    hbnds = implicit_relax_h(h1,hj1,p_mc,pmid,X,P,opts1,param)

    fbnds = implicit_relax_f(f,h1,hj1,X,P,p,pmid,opts1,param)

    fgbnds = implicit_relax_fg(f,g,h1,hj1,X,P,p,pmid,opts1,param)

    param = gen_expansion_params(h1,hj1,X,P,pmid,opts1)

    @test isapprox(param[2][1].cc,-0.4,atol=1E-4)
    @test isapprox(param[2][1].cv,-0.78,atol=1E-4)
    hbnds = implicit_relax_h(h1,hj1,p_mc,pmid,X,P,opts1,param)
    @test isapprox(hbnds[1].cc,-0.5209127114919797,atol=1E-4)
    @test isapprox(hbnds[1].cv,-0.6307841683146562,atol=1E-4)
    @test isapprox(hbnds[1].cc_grad[1],0.0455312,atol=1E-4)
    @test isapprox(hbnds[1].cv_grad[1],0.0941469,atol=1E-4)

    fbnds = implicit_relax_f(f,h1,hj1,X,P,p,pmid,opts1,param)
    @test isapprox(fbnds.cc,3.774523731048122,atol=1E-4)
    @test isapprox(fbnds.cv,2.5572882333553046,atol=1E-4)
    @test isapprox(fbnds.cc_grad[1],0.873187,atol=1E-4)
    @test isapprox(fbnds.cv_grad[1],0.792877,atol=1E-4)

    fgbnds = implicit_relax_fg(f,g,h1,hj1,X,P,p,pmid,opts1,param)
    @test isapprox(fgbnds[1].cc,3.774523731048122,atol=1E-4)
    @test isapprox(fgbnds[1].cv,2.5572882333553046,atol=1E-4)
    @test isapprox(fgbnds[2][1].cc,3.774523731048122,atol=1E-4)
    @test isapprox(fgbnds[2][1].cv,2.5572882333553046,atol=1E-4)
    @test isapprox(fgbnds[2][2].cc,11.274523731048122,atol=1E-4)
    @test isapprox(fgbnds[2][2].cv,10.057288233355305,atol=1E-4)

    @test isapprox(fgbnds[1].cc_grad[1],0.873187,atol=1E-4)
    @test isapprox(fgbnds[1].cv_grad[1],0.792877,atol=1E-4)
    @test isapprox(fgbnds[2][1].cc_grad[1],0.873187,atol=1E-4)
    @test isapprox(fgbnds[2][1].cv_grad[1],0.792877,atol=1E-4)
    @test isapprox(fgbnds[2][2].cc_grad[1],1.87319,atol=1E-4)
    @test isapprox(fgbnds[2][2].cv_grad[1],1.79288,atol=1E-4)
end


@testset "Test Base Implicit Routines (In-place, Single Step)" begin
    EAGO.set_mc_differentiability!(0)
    mc_opts1 = mc_opts(0.5,1,:Dense,:Newton,1,1,1E-10)

    P = [IntervalType(6.0,9.0)]
    X = [IntervalType(-0.78,-0.4)]
    pmid = mid.(P)
    pref_mc = [MC{1}(pmid[1],P[1],1)]
    p_mc = [MC{1}(pmid[1],P[1],1)]

    flt_param = [2.0; 1.0]

    Y = [0.0 0.0; 0.0 0.0]
    H = fill(zero(MC{1}), (1,))
    J = fill(zero(MC{1}), (2,2))

    x_mc = fill(zero(MC{1}), (1,))
    xp_mc = fill(zero(MC{1}), (1,))
    xa_mc = fill(zero(MC{1}), (1,))
    xA_mc = fill(zero(MC{1}), (1,))
    z_mc = fill(zero(MC{1}), (1,))
    aff_mc = fill(zero(MC{1}), (1,))

    interval_bnds = true
    precond = true

    z_mc =[MC{1}(IntervalType(mid(X[1])))]
    t1 = z_mc[1]^2
    t2 = z_mc[1]*p_mc[1]
    t3 = 4.0
    t4 = t1 + t2 + t3

    function h3!(out, z, xp, p, t)
        out[1] = z[1]^2 + z[1]*p[1] + 4.0
    end
    function hj3!(out, z, xp, p, t)
        out[1] = 2.0*z[1]+p[1]
    end

    param = fill(zero(MC{1}), (2,))# probably this

    gen_expansion_params!(h3!, hj3!, pref_mc, xp_mc, x_mc, xa_mc, xA_mc, z_mc, aff_mc,
                          X, P, mc_opts1, param, H, J, Y, interval_bnds, flt_param, precond)
    x_mc_inplace = deepcopy(x_mc)
    params_inplace = deepcopy(param)

    implicit_relax_h!(h3!, hj3!, p_mc, pref_mc, xp_mc, x_mc, xa_mc, xA_mc, z_mc, aff_mc,
                      X, P, mc_opts1, param, H, J, Y, interval_bnds, flt_param, precond)
    hrelax_inplace = deepcopy(x_mc)

    hbnds = x_mc
    @test isapprox(hrelax_inplace[1].cc,-0.5209127114919797,atol=1E-4)
    @test isapprox(hrelax_inplace[1].cv,-0.6307841683146562,atol=1E-4)
    @test isapprox(hrelax_inplace[1].cc_grad[1],0.0455312,atol=1E-4)
    @test isapprox(hrelax_inplace[1].cv_grad[1],0.0941469,atol=1E-4)

    @test isapprox(params_inplace[2][1].cc,-0.4,atol=1E-4)
    @test isapprox(params_inplace[2][1].cv,-0.78,atol=1E-4)
end

@testset "Test Base Implicit Routines (In-place, Two Step)" begin
    EAGO.set_mc_differentiability!(0)
    mc_opts1 = mc_opts(0.5,2,:Dense,:Newton,1,1,1E-10)

    f(x, p) = x[1]*p[1]+p[1]
    g(x, p) = [x[1]*p[1]+p[1];
               x[1]*p[1]+2*p[1]]

    P = [IntervalType(6.0,9.0)]
    X = [IntervalType(-0.78,-0.4)]
    pmid = mid.(P)
    pref_mc = [MC{1}(pmid[1],P[1],1)]
    p_mc = [MC{1}(pmid[1],P[1],1)]

    flt_param = [2.0; 1.0]

    Y = [0.0 0.0; 0.0 0.0]
    H = fill(zero(MC{1}), (1,))
    J = fill(zero(MC{1}), (2,2))

    x_mc = fill(zero(MC{1}), (1,))
    xp_mc = fill(zero(MC{1}), (1,))
    xa_mc = fill(zero(MC{1}), (1,))
    xA_mc = fill(zero(MC{1}), (1,))
    z_mc = fill(zero(MC{1}), (1,))
    aff_mc = fill(zero(MC{1}), (1,))

    interval_bnds = true
    precond = true

    z_mc =[MC{1}(IntervalType(mid(X[1])))]
    t1 = z_mc[1]^2
    t2 = z_mc[1]*p_mc[1]
    t3 = 4.0
    t4 = t1 + t2 + t3

    function h3!(out, z, xp, p, t)
        out[1] = z[1]^2 + z[1]*p[1] + 4.0
    end
    function hj3!(out, z, xp, p, t)
        out[1] = 2.0*z[1]+p[1]
    end

    param = fill(zero(MC{1}), (1,3))

    gen_expansion_params!(h3!, hj3!, pref_mc, xp_mc, x_mc, xa_mc, xA_mc, z_mc, aff_mc,
                          X, P, mc_opts1, param, H, J, Y, interval_bnds, flt_param, precond)
    x_mc_inplace = deepcopy(x_mc)
    params_inplace = deepcopy(param)

    implicit_relax_h!(h3!, hj3!, p_mc, pref_mc, xp_mc, x_mc, xa_mc, xA_mc, z_mc, aff_mc,
                      X, P, mc_opts1, param, H, J, Y, interval_bnds, flt_param, precond)
    hrelax_inplace = deepcopy(x_mc)

    hbnds = x_mc
    @test isapprox(hrelax_inplace[1].cc,-0.5209127114919797,atol=1E-4)
    @test isapprox(hrelax_inplace[1].cv,-0.6307841683146562,atol=1E-4)
    @test isapprox(hrelax_inplace[1].cc_grad[1],0.045531201933640894,atol=1E-4)
    @test isapprox(hrelax_inplace[1].cv_grad[1],0.09414689079323225,atol=1E-4)

    @test isapprox(params_inplace[2][1].cc,-0.4,atol=1E-4)
    @test isapprox(params_inplace[2][1].cv,-0.78,atol=1E-4)
end
