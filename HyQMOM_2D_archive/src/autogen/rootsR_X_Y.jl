"""
    rootsR_X_Y(Y, e11, e12, e13, e23, e24, e34, e44, ex)

Find roots for realizability checking with specific matrix elements.
"""
function rootsR_X_Y(Y, e11, e12, e13, e23, e24, e34, e44, ex)
    @fastmath begin
        t2 = e12^2
        t3 = e13^2
        t4 = e23^2
        t5 = e24^2
        t6 = e34^2
        t7 = ex^2
        t8 = Y*ex*2.0
        t9 = e13*e34*2.0
        t10 = Y*e11*e44
        t11 = 1.0/Y
        t15 = sqrt(3.0)
        t17 = Y*e12*e24*2.0
        t20 = e13*e23*e24*2.0
        t21 = e12*e23*e34*2.0
        t26 = e13*e34*ex*-2.0
        t32 = e11*e23*e24*e34*2.0
        t33 = e12*e13*e23*e44*2.0
        
        t12 = t11^2
        t13 = t11^3
        t14 = -t8
        t16 = Y*t7
        t18 = e11*t6
        t19 = e44*t3
        t23 = -t17
        t24 = -t10
        t25 = ex*t4*2.0
        t27 = Y*e11*t5
        t28 = Y*e44*t2
        t29 = e12*e24*t8
        t30 = e11*e44*t4
        t37 = t3*t5
        t38 = t2*t6
        t39 = t4*t7
        t36 = -t25
        t46 = t4+t9+t14
        t60 = t11*(t27+t28-t29+t30-t32-t33-t37-t38-t39+e12*e13*e24*e34*2.0+e13*e23*e24*ex*2.0+e12*e23*e34*ex*2.0)*(-0.5)
        
        t47 = t46^2
        t48 = t46^3
        t49 = (t11*t46)/3.0
        t53 = t16+t18+t19+t20+t21+t23+t24+t26+t36
        t50 = (t12*t47)/9.0
        t51 = (t13*t48)/27.0
        t54 = (t11*t53)/3.0
        t55 = (t12*t46*t53)/6.0
        t56 = -t55
        t59 = -(t50-t54)^3
        t62 = (-t51+t55+(t11*(t27+t28-t29+t30-t32-t33-t37-t38-t39+e12*e13*e24*e34*2.0+e13*e23*e24*ex*2.0+e12*e23*e34*ex*2.0))/2.0)^2
        
        t61 = t51+t56+t60
        t63 = t59+t62
        t64 = sqrt(t63)
        t65 = t61+t64
        t66 = cbrt(t65)
        t67 = 1.0/t66
        t68 = t66/2.0
        t69 = -t68
        t70 = -t67*(t50-t54)
        t71 = t67*(t50-t54)*(-0.5)
        t72 = t66+t70
        t73 = t15*t72*0.5im
        
        R = [t49+t66+t67*(t50-t54),
             t49+t69+t71-t73,
             t49+t69+t71+t73]
        
        return R
    end
end
