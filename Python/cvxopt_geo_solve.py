import numpy
import time
from cvxopt import solvers, matrix, spmatrix, spdiag, log, exp
class Cvxopt_Solve:
    def __init__(self, marg_xy, marg_xz, _set_to_zero=set()):
        # marg_xy is a dictionary     (x,y) --> positive double
        # marg_xz is a dictionary     (x,z) --> positive double

        self.orig_marg_xy      = None
        self.orig_marg_xz      = None
        self.q_xy              = None
        self.q_xz              = None
        self.var_idx           = None
        self.dual_var_idx      = None
        self.X            = None
        self.Y            = None
        self.Z            = None
        self.F            = None
        self.K            = None
        self.g            = None
        self.marg_of_idx  = None  # triples xyz -- the missing one is None
        self.create_ieqs_called = False
        #self.p_0          = None # initial solution
        self.solver_ret   = None
        self.p_final      = None
        self.set_to_zero  = _set_to_zero
        self.est_opt      = None
        # Actual code:
        self.orig_marg_xy = dict(marg_xy)
        self.orig_marg_xz = dict(marg_xz)
        self.X = set( [ x   for x,y in self.orig_marg_xy.keys() ] + [ x   for x,z in self.orig_marg_xz.keys() ] )
        self.Y = set( [  y  for x,y in self.orig_marg_xy.keys() ] )
        self.Z = set(                                               [  z  for x,z in self.orig_marg_xz.keys() ] )
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

    # create_ieqs():
    def create_ieqs(self):
        # The point is that all the double values are considered non-zero
        # This function
        #  - creates the sets X,Y,Z
        #  - creates the dictionary dual_var_idx: dual variables --> index of the dual variable
        #  - creates the matrix F: 1st column is value two types of marginal equations: p(x,y,*)==p_xy(x,y), p(x,*,z)==p_xz(x,z), the rest is their respective marginal matrix
        #  - creates the matrix g: all zeros
        #  - creates the list K : contains the dimensions of the of F_i's where F^T = [F_0, F_1, ... , F_m]^T
        #  - F has full rank: unneeded ieqs are thrown out
        if self.create_ieqs_called:
            print("Some dork called create_ieqs() twice...")
            exit(1)
        self.create_ieqs_called = True 

        count_dual_vars = 0
        list_K = [] # list of dimesions of F_i's i = |YxZ| + 1
        list_K.append(1) # dimension of F_0 is 1
        self.dual_var_idx    = dict()
        for y in self.Y:
            for z in self.Z:
                count_x = 0 # Dimensions of F_yz : number of x in yz slice
                for x in self.X:
                        if (x,y,z) not in self.set_to_zero  and  (x,y) in self.q_xy.keys() and (x,z) in self.q_xz.keys():
                            self.dual_var_idx[ (x,y,z) ] = count_dual_vars
                            count_dual_vars += 1
                            count_x += 1
                if count_x != 0:
                    list_K.append(count_x)
        tmp_g = matrix(1,(count_dual_vars + 1,1))
        list_b           = [] # list of RHSs
        list_Ft          = [] # coefficient matrix in row-major (need column-major later)
        list_Ft_throwout = [] # DEBUG: omited equations go here---then we check whether the ranks are equal
        numo_thrownout   = 0

        self.marg_of_idx = []
        count_vars = 0
        self.var_idx = dict()

        # Create variables' dictionary from marginals 

        # xy-marginal equations
        for xy,rhs in self.q_xy.items():
            x,y = xy
            f = [ 0   for xyz in self.dual_var_idx.keys() ] # initialize the whole row with 0
            for z in self.Z:
                if (x,y,z) in self.dual_var_idx.keys():
                    i  = self.dual_var_idx[ (x,y,z) ]
                    f[i] = 1.
            # We test if adding this equation increases the rank.
            # Because of the deleted variables (in set_to_zero), I don't know of a better way to do this...
            tmp_Ft = matrix( list_Ft+f, ( len(self.dual_var_idx),  len(list_b)+1 ), 'd' )
            if numpy.linalg.matrix_rank( tmp_Ft ) > len(list_b):
                list_Ft += f # splice !
                list_b.append( rhs )
                self.marg_of_idx.append(  (x,y,None)  )
                self.var_idx[(x,y,None)] = count_vars
                count_vars += 1

        # xz-marginal equations
        for xz,rhs in self.q_xz.items():
            x,z = xz
            f = [ 0   for xyz in self.dual_var_idx.keys() ] # initialize the whole row with 0
            for y in self.Y:
                if (x,y,z) in self.dual_var_idx.keys():
                    i = self.dual_var_idx[ (x,y,z) ]
                    f[i] = 1.
            # Rank-check again:
            tmp_Ft = matrix( list_Ft+f, ( len(self.dual_var_idx),  len(list_b)+1 ), 'd' )
            if numpy.linalg.matrix_rank( tmp_Ft ) > len(list_b):
                list_Ft += f # splice !
                list_b.append( rhs )
                self.marg_of_idx.append(  (x,None,z)  )
                self.var_idx[(x,None,z)] = count_vars
                count_vars += 1
        # Now we create the CvxOpt gp matrix.
        Ft      = matrix( list_Ft,                  ( len(self.dual_var_idx),  len(list_b)                ), 'd' )
        tmp_F   = Ft.T
        rk = numpy.linalg.matrix_rank(tmp_F)
        if ( rk != len(list_b) ):
            print("BUG: There's something wrong with the rank of the coefficient matrix: it is ",rk," it should be ",len(list_b))
            exit(1)
        tmp_b = matrix(list_b, (len(list_b), 1), 'd')
        self.F = matrix([[tmp_b],[tmp_F]]).T
        # The partioning dimesnions of the Cvxopt matrix + reduntant scalars g. 
        self.K = list_K
        self.g = log(tmp_g)
        # Print if neccessarly 
        print("F: ",self.F)
        print("K: ",self.K)
        print("g: ",self.g)
        # dim_space = len(self.var_idx)-rk
        # print("Solution space has dimension ",dim_space)
    #^ create_ieqs()
    
    # solve_it():
    def solve_it(self):
        self.q_xy = self.tidy_up_distrib(self.orig_marg_xy)
        self.q_xz = self.tidy_up_distrib(self.orig_marg_xz)

        self.create_ieqs()
        start_opt = time.clock()
        #solvers.options['maxiters']= 200
        self.solver_ret   = solvers.gp(K=self.K, F=self.F, g=self.g)
        print("Solver terminated with status ",self.solver_ret['status'])
        self.est_opt = (time.clock() - start_opt)
    #^ solve_it()

    def callback(self,p=None, zz=None):
        N = len(self.dual_var_idx)
        if p is None:
            list_p_0 = [ 0.   for xyz in self.dual_var_idx.keys() ]
            for xyz,i in self.dual_var_idx.items():
                list_p_0[i] = 1. # self.p_0[xyz]
                # This is returns the starting solution for the iterative solution of the CP --- this is the 1st point to experiment with other distributions with the same marginals.
                return 0, matrix(list_p_0, (N,1), 'd' )

        # check if p is in the feasible region for the objective function:
        if min(p) <= 0 or max(p) > 1:
            return None

        p_dict = dict( (xyz,p[i]) for xyz,i in self.dual_var_idx.items() )
        p_yz = marginal_yz(p_dict)

        # Compute f(p)
        f = 0
        for xyz,i in self.dual_var_idx.items():
            x,y,z = xyz
            if p[i] > 0: f += p[i]*log(p[i]/p_yz[y,z])
        return f
    #^ callback()

    def check_feasibility(self,pdf):

        # Get the number of variables 
        var_num = len(self.var_idx)

        # Get the number of dual variables
        dual_var_num = len(self.dual_var_idx)
        x_sz = len(self.X)
        y_sz = len(self.Y)
        z_sz = len(self.Z)

        # Compute marginals 
        p_xy = marginal_xy(pdf)
        p_xz = marginal_xz(pdf)
        print(p_xy)
        # Compute |Y_x|
        c_xy = 0
        for x in self.X:
            for y in self.Y:
                if (x,y) in p_xy.keys() and p_xy[(x,y)] > 0:
                    c_xy += 1

        # Compute |Z_x|
        c_xz = 0
        for x in self.X:
            for z in self.Z:
                if (x,z) in p_xz.keys() and p_xy[(x,z)] > 0:
                    c_xz += 1
        dual_feas_dim = dual_var_num + x_sz - c_xy - c_xz
        status = self.solver_ret['status'] # This is a string

        # Make sol_list
        print(self.solver_ret['x'])
        print(self.var_idx)
        sol = numpy.zeros((len(self.var_idx),1))
        iter = 0
        for i in self.var_idx.values():
            sol[iter] = self.solver_ret['x'][i]
            iter += 1
        
        sol_list = []
        for i in self.var_idx.values():
            sol_list.append(self.solver_ret['x'][i])

        # Make dual_sol ( q prob distribution ) 
        q = numpy.zeros((len(self.dual_var_idx),1))
        iter = 0
        for i in self.dual_var_idx.values():
            q[iter] = self.solver_ret['znl'][i]
            iter += 1
        
        q_list = []
        for i in self.dual_var_idx.values():
            q_list.append(self.solver_ret['znl'][i])
        
        dual_obj_val = self.callback(q_list)

        q_nonneg_viol = max(-min(q),0)
        q_min_entry   = max(min(q),0)

        print("q: ", q)
        print("Dual obj val: ", dual_obj_val)
        print("q_noneg_viol: ", q_nonneg_viol)
        print("q_min_entry: ", q_min_entry)
        
        # self.A*p - self.b

        # equation = numpy.matmul(q.T,self.A.T).T - self.b

        # marginals_1   = numpy.linalg.norm(equation, 1)
        # marginals_2   = numpy.linalg.norm(equation, 2)
        # marginals_Inf = numpy.linalg.norm(equation, numpy.inf)

        # llambda = self.solver_ret['y'] # array; fits with A I guess...

        # mu = gradient.T + numpy.matmul(self.A.T,llambda)
        # # mu_nonneg_viol
        # mu_nonneg_viol = -min(mu)
        
        # complementarity_max = max( numpy.multiply( numpy.absolute(mu), numpy.absolute(q) ) )
        # complementarity_sum = sum( numpy.multiply( numpy.absolute(mu), numpy.absolute(q) ) )

        CI   = -1.0
        SI   = -1.0
        UI_Y = -1.0
        UI_Z = -1.0

        opt_time        = self.est_opt

        return  var_num, dual_var_num, x_sz, y_sz, z_sz, dual_feas_dim, status, dual_obj_val, q_nonneg_viol, q_min_entry[0], CI, SI, UI_Y, UI_Z, opt_time
    #^ check_feasibility()
    
    def do_it(self, pdf):

        self.solve_it()

        var_num, dual_var_num, x_sz, y_sz, z_sz, dual_feas_dim, status, dual_obj_val, q_nonneg_viol, q_min_entry[0], CI, SI, UI_Y, UI_Z, opt_time = self.check_feasibility(pdf)
        
        return  var_num, dual_var_num, x_sz, y_sz, z_sz, dual_feas_dim, status, dual_obj_val, q_nonneg_viol, q_min_entry[0], CI, SI, UI_Y, UI_Z, opt_time
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

    var_num, dual_var_num, x_sz, y_sz, z_sz, dual_feas_dim, status, dual_obj_val, q_nonneg_viol, q_min_entry[0], CI, SI, UI_Y, UI_Z, opt_time = cvx.do_it(pdf)
    
    return  var_num, dual_var_num, x_sz, y_sz, z_sz, dual_feas_dim, status, dual_obj_val, q_nonneg_viol, q_min_entry[0], CI, SI, UI_Y, UI_Z, opt_time

###########
# Test Run
###########
pdf = {(0,0,0):.5, (1,1,1):.5}
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

solve_PDF(pdf)
