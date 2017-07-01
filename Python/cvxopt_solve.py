import numpy
import time
from cvxopt import solvers, matrix, spmatrix, spdiag, log
class Cvxopt_Solve:
    def __init__(self, marg_xy, marg_xz, _set_to_zero=set()):
        # marg_xy is a dictionary     (x,y) --> positive double
        # marg_xz is a dictionary     (x,z) --> positive double

        self.orig_marg_xy = None
        self.orig_marg_xz = None
        self.q_xy         = None
        self.q_xz         = None
        self.var_idx      = None
        self.X            = None
        self.Y            = None
        self.Z            = None
        self.A            = None
        self.b            = None
        self.marg_of_idx  = None  # triples xyz -- the missing one is None
        self.create_equations_called = False
        self.G            = None
        self.h            = None
        self.create_ieqs_called = False
        self.p_0          = None # initial solution
        self.solver_ret   = None
        self.p_final      = None
        self.set_to_zero  = _set_to_zero
        self.est_opt      = None
        self.num_eval_f      = None
        self.num_eval_grad_f = None
        self.num_eval_g      = None
        self.num_eval_jac_g  = None
        self.num_hessevals   = None
        # Actual code:
        self.orig_marg_xy = dict(marg_xy)
        self.orig_marg_xz = dict(marg_xz)
        self.X = set( [ x   for x,y in self.orig_marg_xy.keys() ] + [ x   for x,z in self.orig_marg_xz.keys() ] )
        self.Y = set( [  y  for x,y in self.orig_marg_xy.keys() ] )
        self.Z = set(                                               [  z  for x,z in self.orig_marg_xz.keys() ] )
        self.num_eval_f      = 0
        self.num_eval_grad_f = 0
        self.num_eval_g      = -1
        self.num_eval_jac_g  = -1
        self.num_hessevals   = 0
        self.est_opt         = 0
    # __init__()

    # tidy_up_distrib():
    def tidy_up_distrib(self,p):
        eps = 0
        # returns a tidied-up copy of p: every entry smaller than eps is treated as 0.
        p_new = dict()
        one = 0.
        for x,r in p.items():
            if r>=eps:
                p_new[x] = r
                one += r
        # Re-normalize --- Should I drop this ??? ?!?
        for x,r in p_new.items():
            p_new[x] = r / one;

        return p_new
    #^ tidy_up_distrib()

    # create_equations():
    def create_equations(self):
        # The point is that all the double values are considered non-zero
        # This function
        #  - creates the sets X,Y,Z
        #  - creates the dictionary var_idx: variables --> index of the variable
        #  - creates the matrices A,b: two types of marginal equations: p(x,y,*)==p_xy(x,y), p(x,*,z)==p_xz(x,z)
        #  - A has full rank: unneeded eqns are thrown out
        if self.create_equations_called:
            print("Some dork called create_equations twice...")
            exit(1)
        self.create_equations_called = True

        count_vars = 0
        self.var_idx    = dict()
        for x in self.X:
            for y in self.Y:
                if (x,y) in self.q_xy.keys():
                    for z in self.Z:
                        if (x,y,z) not in self.set_to_zero  and  (x,z) in self.q_xz.keys():
                            self.var_idx[ (x,y,z) ] = count_vars
                            count_vars += 1

        list_b           = [] # list of RHSs
        list_At          = [] # coefficient matrix in row-major (need column-major later)
        list_At_throwout = [] # DEBUG: omited equations go here---then we check whether the ranks are equal
        numo_thrownout   = 0

        self.marg_of_idx = []
        # xy-marginal equations
        for xy,rhs in self.q_xy.items():
            x,y = xy
            a = [ 0   for xyz in self.var_idx.keys() ] # initialize the whole row with 0
            for z in self.Z:
                if (x,y,z) in self.var_idx.keys():
                    i  = self.var_idx[ (x,y,z) ]
                    a[i] = 1.
            # list_At += a # splice !
            # list_b.append( rhs )
            # self.marg_of_idx.append(  (x,y,None)  )

            # We test if adding this equation increases the rank.
            # Because of the deleted variables (in set_to_zero), I don't know of a better way to do this...
            tmp_At = matrix( list_At+a, ( len(self.var_idx),  len(list_b)+1 ), 'd' )
            if numpy.linalg.matrix_rank( tmp_At ) > len(list_b):
                list_At += a # splice !
                list_b.append( rhs )
                self.marg_of_idx.append(  (x,None,z)  )

        # xz-marginal equations
        for xz,rhs in self.q_xz.items():
            x,z = xz
            a = [ 0   for xyz in self.var_idx.keys() ] # initialize the whole row with 0
            for y in self.Y:
                if (x,y,z) in self.var_idx.keys():
                    i = self.var_idx[ (x,y,z) ]
                    a[i] = 1.
            # Rank-check again:
            tmp_At = matrix( list_At+a, ( len(self.var_idx),  len(list_b)+1 ), 'd' )
            if numpy.linalg.matrix_rank( tmp_At ) > len(list_b):
                list_At += a # splice !
                list_b.append( rhs )
                self.marg_of_idx.append(  (x,None,z)  )

        # Now we create the CvxOpt matrix.
        self.b  = matrix( list_b,                   ( len(list_b),        1                          ), 'd' )
        At      = matrix( list_At,                  ( len(self.var_idx),  len(list_b)                ), 'd' )
        self.A = At.T
        rk = numpy.linalg.matrix_rank(self.A)
        if ( rk != len(list_b) ):
            print("BUG: There's something wrong with the rank of the coefficient matrix: it is ",rk," it should be ",len(list_b))
            exit(1)
        dim_space = len(self.var_idx)-rk
        print("Solution space has dimension ",dim_space)
    #^ create_equations()

    # create_ieqs():
    def create_ieqs(self):
        if not self.create_equations_called:
            print("You have to call create_equations() before calling create_ieqs()")
            exit(1)
        if self.create_ieqs_called:
            print("Some dork called create_ieqs() twice...")
            exit(1)
        self.create_ieqs_called = True
        self.G = spdiag( matrix( -1., (len(self.var_idx),1), 'd' ) )
        self.h = matrix( 0, (len(self.var_idx),1), 'd' )
    #^ create_ieqs()

    # CALLBACK fn for computing  f, grad f, Hess f
    def callback(self,p=None, zz=None):
        N = len(self.var_idx)
        if p is None:
            list_p_0 = [ 0.   for xyz in self.var_idx.keys() ]
            for xyz,i in self.var_idx.items():
                list_p_0[i] = 1. # self.p_0[xyz]
            # This is returns the starting solution for the iterative solution of the CP --- this is the 1st point to experiment with other distributions with the same marginals.
            return 0, matrix(list_p_0, (N,1), 'd' )

        # check if p is in the feasible region for the objective function:
        if min(p) <= 0 or max(p) > 1:
            return None

        p_dict = dict( (xyz,p[i]) for xyz,i in self.var_idx.items() )
        p_yz = marginal_yz(p_dict)

        # Compute f(p)
        f = 0
        for xyz,i in self.var_idx.items():
            x,y,z = xyz
            if p[i] > 0: f += p[i]*log(p[i]/p_yz[y,z])

        # Compute gradient-transpose Df(p)
        list_Df = [ 0. for xyz in self.var_idx.keys() ]
        for xyz,i in self.var_idx.items():
            x,y,z = xyz
            pyzyz = p_yz[y,z]
            if p[i] > 0:     list_Df[i] = log( p[i] / pyzyz )
            elif pyzyz <= 0: list_Df[i] = -log(len(self.X))
            else:            list_Df[i] = -exp(max(10,len(self.X)))
        Df = matrix(list_Df, (1,N), 'd')

        if zz is None:
            self.num_eval_f      += 1
            self.num_eval_grad_f += 1
            return f,Df

        # Compute zz[0] * Hess f
        # This will be a sparse matrix
        entries = []
        rows    = []
        columns = []
        for xyz,i in self.var_idx.items():
            x,y,z = xyz
            p_yz__x = p_yz[y,z] - p[i] # sum_{* \ne x} p(*,y,z).
            for x_ in self.X:
                if x_==x: # diagonal
                    rows.append( i )
                    columns.append( i )
                    tmp_quot = zz[0] * p_yz__x / p_yz[y,z] # 1/p[x,y,z] - 1/p[*,y,z] = ( p[*,y,z] - p[x,y,z] )/( p[*,y,z] p[x,y,z] )
                    if p[i] > 0:   entries.append( tmp_quot / p[i] )
                    else:
                        print("TROUBLE computing Hessian (diagonal)")
                        entries.append( -1.e-300 )
                else: # off diagonal
                    if (x_,y,z) in self.var_idx:
                        j = self.var_idx[ (x_,y,z) ]
                        val = - zz[0] / p_yz[y,z]
                        rows.append( i )
                        columns.append( j )
                        entries.append( val )
                # if diagonal
            # for x_
        # for xyz,i

        zH = spmatrix( entries, rows, columns, (N,N), 'd')
        # if self.verbose_output: print("p=",list(p))
        self.num_eval_f      += 1
        self.num_eval_grad_f += 1
        self.num_hessevals   += 1
        return f,Df,zH
    #^ callback()


    # make_initial_solution()
    def make_initial_solution(self):
        self.p_0 = dict()
        for xyz in self.var_idx.keys():
            self.p_0[xyz] = 1.
    #^ make_initial_solution()

    # solve_it():
    def solve_it(self):
        self.q_xy = self.tidy_up_distrib(self.orig_marg_xy)
        self.q_xz = self.tidy_up_distrib(self.orig_marg_xz)

        self.create_equations()
        self.create_ieqs()
        self.make_initial_solution()
        start_opt = time.clock()
        self.solver_ret   = solvers.cp(self.callback, G=self.G, h=self.h, A=self.A, b=self.b)
        print("Solver terminated with status ",self.solver_ret['status'])
        self.est_opt = (time.clock() - start_opt)
    #^ solve_it()

    def check_feasibility(self):
        var_num = len(self.var_idx)
        x_sz = len(self.X)
        y_sz = len(self.Y)
        z_sz = len(self.Z)
        

        status = self.solver_ret['status'] # This is a string

        # Make q
        q = numpy.zeros((len(self.var_idx),1))
        iter = 0
        for i in self.var_idx.values():
            q[iter] = self.solver_ret['x'][i]
            iter += 1
        
        # q_ = dict()
        # for xyz,i in self.var_idx.items():
        #     q_[xyz] = self.solver_ret['x'][i]
        # self.p_final = q_
        q_list = []
        for i in self.var_idx.values():
            q_list.append(self.solver_ret['x'][i])
        
        obj_val,gradient = self.callback(q_list)

        q_nonneg_viol = max(-min(q),0)
        q_min_entry   = max(min(q),0)
        
        # self.A*p - self.b

        equation = numpy.matmul(q.T,self.A.T).T - self.b

        marginals_1   = numpy.linalg.norm(equation, 1)
        marginals_2   = numpy.linalg.norm(equation, 2)
        marginals_Inf = numpy.linalg.norm(equation, numpy.inf)

        llambda = self.solver_ret['y'] # array; fits with A I guess...
        # print(numpy.matmul(self.A.T,llambda))
        # print(gradient.T)
        # print("lambda")
        # print(llambda)
        mu = gradient.T + numpy.matmul(self.A.T,llambda)
        # mu_nonneg_viol
        mu_nonneg_viol = -min(mu)
        
        complementarity_max = max( numpy.multiply( numpy.absolute(mu), numpy.absolute(q) ) )
        complementarity_sum = sum( numpy.multiply( numpy.absolute(mu), numpy.absolute(q) ) )

        CI   = -1.0
        SI   = -1.0
        UI_Y = -1.0
        UI_Z = -1.0

        num_eval_f      = self.num_eval_f
        num_eval_grad_f = self.num_eval_grad_f
        num_eval_g      = self.num_eval_g
        num_eval_jac_g  = self.num_eval_jac_g
        num_hessevals   = self.num_hessevals
        opt_time        = self.est_opt

        return  var_num, x_sz, y_sz, z_sz, status, obj_val, q_nonneg_viol, q_min_entry[0], marginals_1, marginals_2, marginals_Inf, mu_nonneg_viol[0], complementarity_max[0], complementarity_sum[0], CI, SI, UI_Y, UI_Z, num_eval_f, num_eval_grad_f, num_eval_g, num_eval_jac_g, num_hessevals, opt_time
    #^ check_feasibility()
    
    def do_it(self):

        self.solve_it()

        var_num, x_sz, y_sz, z_sz, status, obj_val, q_nonneg_viol, q_min_entry, marginals_1, marginals_2, marginals_Inf, mu_nonneg_viol, complementarity_max, complementarity_sum, CI, SI, UI_Y, UI_Z, num_eval_f, num_eval_grad_f, num_eval_g, num_eval_jac_g, num_hessevals, opt_time = self.check_feasibility()
        
        return  var_num, x_sz, y_sz, z_sz, status, obj_val, q_nonneg_viol, q_min_entry, marginals_1, marginals_2, marginals_Inf, mu_nonneg_viol, complementarity_max, complementarity_sum, CI, SI, UI_Y, UI_Z, num_eval_f, num_eval_grad_f, num_eval_g, num_eval_jac_g, num_hessevals, opt_time
    #^do_it()

    
#^ class Cvxopt_Solve

###############
# USE IT
###############


# Marginals

def marginal_xy(p):
    marg = dict()
    for xyz,r in p.items():
        x,y,z = xyz
        if (x,y) in marg.keys():    marg[(x,y)] += r
        else:                       marg[(x,y)] =  r
    return marg

def marginal_xz(p):
    marg = dict()
    for xyz,r in p.items():
        x,y,z = xyz
        if (x,z) in marg.keys():   marg[(x,z)] += r
        else:                      marg[(x,z)] =  r
    return marg


def marginal_yz(p):
    marg = dict()
    for xyz,r in p.items():
        x,y,z = xyz
        if (y,z) in marg.keys():    marg[(y,z)] += r
        else:                       marg[(y,z)] =  r
    return marg

def solve_PDF(pdf):
    p_xy = marginal_xy(pdf)
    p_xz = marginal_xz(pdf)
    cvx = Cvxopt_Solve(p_xy,p_xz)
    var_num, x_sz, y_sz, z_sz, status, obj_val, q_nonneg_viol, q_min_entry, marginals_1, marginals_2, marginals_Inf, mu_nonneg_viol, complementarity_max, complementarity_sum, CI, SI, UI_Y, UI_Z, num_eval_f, num_eval_grad_f, num_eval_g, num_eval_jac_g, num_hessevals, opt_time = cvx.do_it()
    
    return var_num, x_sz, y_sz, z_sz, status, obj_val, q_nonneg_viol, q_min_entry, marginals_1, marginals_2, marginals_Inf, mu_nonneg_viol, complementarity_max, complementarity_sum, CI, SI, UI_Y, UI_Z, num_eval_f, num_eval_grad_f, num_eval_g, num_eval_jac_g, num_hessevals, opt_time


###########
# Test Run
###########
#pdf = {(0,0,0):.5, (1,1,1):.5}
# pdf = {((0,0),0,0):.25, ((0,1),0,1):.25, ((1,0),1,0):.25, ((1,1),1,1):.25}
# pdf = {(0,0,0):.25,(1,0,1):.25,(1,1,0):.25,(0,1,1):.25}
# pdf = {(0,0,0):.25, (0,0,1):.25, (0,1,0):.25, (1,1,1):.25}
# pdf = {
#             ((0,0), (0,0), (0,0)): .125,
#             ((1,0), (0,0), (1,0)): .125,
#             ((1,0), (1,0), (0,0)): .125,
#             ((0,0), (1,0), (1,0)): .125,
#             ((0,1), (0,1), (0,1)): .125,
#             ((1,1), (0,1), (1,1)): .125,
#             ((1,1), (1,1), (0,1)): .125,
#             ((0,1), (1,1), (1,1)): .125
#         }
# pdf = {
#             ((0,0,0,0), (0,0), (0,0)): 1/32,
#             ((1,0,0,0), (0,0), (1,0)): 1/32,
#             ((0,0,1,0), (0,0), (2,0)): 1/32,
#             ((1,0,1,0), (0,0), (3,0)): 1/32,
#             ((1,0,0,0), (1,0), (0,0)): 1/32,
#             ((0,0,0,0), (1,0), (1,0)): 1/32,
#             ((1,0,1,0), (1,0), (2,0)): 1/32,
#             ((0,0,1,0), (1,0), (3,0)): 1/32,
#             ((0,1,0,0), (2,0), (0,0)): 1/32,
#             ((1,1,0,0), (2,0), (1,0)): 1/32,
#             ((0,1,1,0), (2,0), (2,0)): 1/32,
#             ((1,1,1,0), (2,0), (3,0)): 1/32,
#             ((1,1,0,0), (3,0), (0,0)): 1/32,
#             ((0,1,0,0), (3,0), (1,0)): 1/32,
#             ((1,1,1,0), (3,0), (2,0)): 1/32,
#             ((0,1,1,0), (3,0), (3,0)): 1/32,
#             ((0,0,0,1), (0,1), (0,1)): 1/32,
#             ((1,0,0,1), (0,1), (1,1)): 1/32,
#             ((0,0,1,1), (0,1), (2,1)): 1/32,
#             ((1,0,1,1), (0,1), (3,1)): 1/32,
#             ((1,0,0,1), (1,1), (0,1)): 1/32,
#             ((0,0,0,1), (1,1), (1,1)): 1/32,
#             ((1,0,1,1), (1,1), (2,1)): 1/32,
#             ((0,0,1,1), (1,1), (3,1)): 1/32,
#             ((0,1,0,1), (2,1), (0,1)): 1/32,
#             ((1,1,0,1), (2,1), (1,1)): 1/32,
#             ((0,1,1,1), (2,1), (2,1)): 1/32,
#             ((1,1,1,1), (2,1), (3,1)): 1/32,
#             ((1,1,0,1), (3,1), (0,1)): 1/32,
#             ((0,1,0,1), (3,1), (1,1)): 1/32,
#             ((1,1,1,1), (3,1), (2,1)): 1/32,
#             ((0,1,1,1), (3,1), (3,1)): 1/32
#         }

# pdf = {
#             ((0,0), 0, 0): .25,
#             ((1,0), 0, 1): .25,
#             ((1,0), 1, 0): .25,
#             ((0,1), 1, 1): .25,
#         }

#solve_PDF(pdf)
