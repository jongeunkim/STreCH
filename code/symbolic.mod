
##### Options
option randseed 1 ;
#option substout 1 ;
option solver scipampl ;

##### model parameters
param constr_vdef_test_leaf default 0 ;


##### Parameters for problem information
param	func symbolic ;

##### Paramteres for tree structure 
set 	NODE ordered	default 1..15 ;
param 	D				:= floor(log(last(NODE))/log(2)) ;
set 	NLEAF 			:= {n in NODE: {2*n, 2*n+1} within NODE} ;
set 	LEAF  			:= {n in NODE: n not in NLEAF} ; 
#display D, NODE, NLEAF, LEAF ; 

##### Paramters for operators and operands
set 	OPER_AVAIL_B	:= {"PLUS", "MINUS", "MULT", "DIV"} ;
set 	OPER_AVAIL_U	:= {"SQR", "SQRT", "EXP", "LOG"} ;
set 	OPER_AVAIL		:= OPER_AVAIL_B union OPER_AVAIL_U ;

set		OPER_NOT_IN_USE within OPER_AVAIL		default { } ;	# set of operators not to include
set 	OPER_IN_USE		within OPER_AVAIL		default OPER_AVAIL ;

set 	OPER_B			:= ( OPER_IN_USE diff OPER_NOT_IN_USE ) inter OPER_AVAIL_B ;
set 	OPER_U			:= ( OPER_IN_USE diff OPER_NOT_IN_USE ) inter OPER_AVAIL_U ;
set 	OPER_BU			:= OPER_B union OPER_U ;

set		VAR ordered		default { } ;	# independent variables "X
set		VAR_D1EXP		default { } ;

set		OPER_V			:= VAR union VAR_D1EXP ;
set		OPER_L			:= if constr_vdef_test_leaf == 0 then OPER_V union {'CST'} else {'LEAF'} ;
set		OPER			:= OPER_BU union OPER_L ;

#display OPER, OPER_B, OPER_U, OPER_L, OPER_V ;


##### Parameters for samples
param	nOBS > 0 integer 							default 10 ;				# the number of observations
set 	OBS ordered 								default 1..nOBS ;
param	loss_weight {OBS}							default 1 ;
param 	x {OBS, VAR union VAR_D1EXP union {'Z'}} 	default 0 ;
param	x_lo {VAR}									default -1e+09 ;
param	x_up {VAR}									default 1e+09 ;
param	noise_level integer							default 0 ;

##### Parameters for defining the feasible space
param	NUM_ACTIVE	integer		default -1 ;
param	NUM_OPER 	integer		default -1 ; 
param	NUM_CST 	integer 	default -1 ;

param 	v_up 		default 100 ;
param 	v_lo 		default -100 ;
param 	cst_up 		default 10 ;
param 	cst_lo 		default -10 ;
param	coeff_up	default 1 ;
param	coeff_lo 	default -1 ;
param	mult_eps	default 0.1 ;
param 	div_eps 	default 0.1 ;
param 	sqrt_eps	default 0.1 ;
param	cst_eps		default 0.1 ;
#param	exp_scale := v_up / exp(v_up);
#param	log_scale := v_up / log(v_up);
param	exp_scale 	default 1;
param	log_scale 	default 1;
param 	log_eps 	:= exp(v_lo / log_scale);
param 	x_eps 		:= 0;

##### Parameters for result display (result.ampl) #####
param inactive integer ;
param EPS_FEAS := 1e-09 ;
param gap ;
param max_feas_gap ;
param param1 ;
param param2 ;
param sym1 symbolic ;


##### DECISON VARIABLE #####
set YSET := {n in NODE, o in OPER: n in NLEAF or o in OPER_L} ;
param yfix {YSET} default -1 ; # -1: nothing, 0 or 1: fix the value to be 0 or 1
param noper_lo 	{o in OPER} 	default 0 ;
param noper_up 	{o in OPER} 	default card(NODE) ; 
param noper_val {o in OPER} 	default -1 ;

var y {YSET}	binary;
var v {i in OBS, n in NODE} >= v_lo, <= v_up ;
var c {n in NODE} >= cst_lo, <= cst_up ;
var p {n in NODE, o in OPER_V} >= coeff_lo, <= coeff_up ;
var e {i in OBS} ;
var noper {o in OPER} >= noper_lo[o], <= noper_up[o] ;

#param ysol {NODE, OPER} default -1 ;

##### OBJECTIVE: minimize MSE #####
param objlambda 	default 0 ; 
minimize Loss: 
	( sum {i in OBS} loss_weight[i] * e[i]^2 ) 
	+ objlambda  * ( sum {(n,o) in YSET}   y[n,o] )
;

subject to errdef {i in OBS}:
	e[i] == v[i,1] - x[i,'Z'] ;

##### Fix solution #####
#subject to Fixysol {n in NODE, o in OPER: yfix[n,o] == 0 or yfix[n,o] == 1}:
#	y[n,o] == yfix[n,o] ;  


##### CONSTR: CONTROL #####

subject to LimitNumCst:
	sum {(n,o) in YSET: o == 'CST'} y[n,o] <= (if NUM_CST > -1 then NUM_CST else 2^(D-1)) ;

subject to LimitNumOper:
	sum {n in NLEAF, o in OPER_BU} y[n,o] <= (if NUM_OPER > -1 then NUM_OPER else card(NLEAF)) ;

subject to LimitNumActive:
	sum {(n,o) in YSET} y[n,o] <= (if NUM_ACTIVE > -1 then NUM_ACTIVE else card(NODE)) ;

##### CONSTR: def noper #####
subject to DefNOper {o in OPER}:
	sum {(n,oo) in YSET: oo == o} y[n,o] == noper[o] ;



##### CONSTR: Linear constraints only with y #####
# Note that if a unary operator is active, then the right child is active.
# It is opposite from [Cozad] and [Neuman].
# This enables to introduce tighter redundancy elimination constraints later.
#subject to FixYRoot {o in OPER_U union OPER_L}:				y[1,o] == 0 ;
#subject to FixYLeaf {n in LEAF, o in OPER_BU}:				y[n,o] == 0 ;
#subject to FixYNotUse {n in NODE, o in OPER_NOT_IN_USE }:		y[n,o] == 0 ;
subject to FixManual 	{(n,o) in YSET: yfix[n,o] in {0,1}}:	y[n,o] == yfix[n,o] ;
subject to FixD1EXP		{(n,o) in YSET: n in LEAF and o in VAR_D1EXP}:	y[n,o] == 0 ;


#subject to Test1 {n in NLEAF, o in {'MULT', 'DIV'}}:			y[n,o] == 0.5 ; 
#subject to Test2 {n in LEAF }:						y[n,'CST'] == 0.5 ;


#subject to Fixyzero {n in NODE, o in OPER: (n in LEAF and o in OPER_BU) or (o in OPER_NOT_IN_USE) }:
#	y[n, o] == 0 ;

# [Cozad]-(1d) at most one operator is active
subject to ActiveAtMost {n in NODE}:
	sum {(nn,o) in YSET: nn == n} y[n,o] <= 1;
	
# [Cozad]-(1i) at least one independent variable is active
subject to ActiveX:
	sum{(n,o) in YSET: o in OPER_V union {'LEAF'}} y[n,o] >= 1;

param constr_ydef default 1 ;
### Option 1: Cozad's formulation
# [Cozad]-(1e) if binary/unary is active, then rch is active
subject to ActiveRch {n in NLEAF: constr_ydef == 1}:
	sum {o in OPER_BU} y[n,o] <= sum {o in OPER} y[(2*n+1),o];
	
# [Cozad]-(1f) if binary is active, then lch is active
subject to ActiveLch {n in NLEAF: constr_ydef == 1}:
	sum {o in OPER_B} y[n,o] <= sum {o in OPER} y[2*n,o];

# [Cozad]-(1h)
subject to InactiveRch {n in NLEAF: constr_ydef == 1}:
	sum {o in OPER_L} y[n,o] <= 1 - sum {o in OPER} y[2*n+1,o];

# [Cozad]-(1g)
subject to InactiveLch {n in NLEAF: constr_ydef == 1}:
	sum {o in OPER_L union OPER_U} y[n,o] <= 1 - sum {o in OPER} y[2*n,o];
	
### [NEW] Option 2: improvement of (1e-1h). Also, improvement of [Newman]'s (18-23)
subject to ActiveInactiveRch {nn in NLEAF: constr_ydef == 2}:
	sum {(n,o) in YSET: n == nn and o in OPER_BU} y[n,o] == sum {(n,o) in YSET: n==2*nn+1} y[n,o];

subject to ActiveInactiveLch {nn in NLEAF: constr_ydef == 2}:
	sum {(n,o) in YSET: n == nn and o in OPER_B}  y[n,o] == sum {(n,o) in YSET: n==2*nn} y[n,o]; 


##### CONSTR: Redundancy elimination #####
param constr_redundancy default 1;
#param constr_vdef_test_leaf default 0 ;

# [Cozad] (2a) no constant binary operation
subject to Red_CstCst {n in NLEAF: constr_redundancy >= 1 and constr_vdef_test_leaf == 0}:
	y[2*n,'CST'] + y[(2*n+1),'CST'] <= 1;

# [Cozad] (2b) no unary of constant
subject to Red_UnaryCst {n in NLEAF: constr_redundancy == 1 and constr_vdef_test_leaf == 0}:
	( sum {o in OPER_U} y[n,o] ) + y[(2*n+1),'CST'] <= 1;

# [Cozad] (3a) no minus constant
subject to Red_MinusCst {n in NLEAF: constr_redundancy == 1 and constr_vdef_test_leaf == 0}:
	y[n,'MINUS']+ y[(2*n+1),'CST'] <= 1;

# [Cozad] (3b) no division by constant
subject to Red_DivCst {n in NLEAF: constr_redundancy == 1 and constr_vdef_test_leaf == 0}:
	y[n,'DIV'] + y[(2*n+1),'CST'] <= 1;

# [Cozad] Nested operations
subject to Red_SqrtSqr {n in NLEAF: constr_redundancy >= 1 and {'SQR', 'SQRT'} within OPER }:
	y[n, 'SQRT'] + y[(2*n+1), 'SQR'] <= 1;

subject to Red_SqrSqrt {n in NLEAF: constr_redundancy >= 1 and {'SQR', 'SQRT'} within OPER }:
	y[n, 'SQR'] + y[(2*n+1), 'SQRT'] <= 1;

subject to Red_ExpLog {n in NLEAF: constr_redundancy >= 1 and 2*n+1 in NLEAF and {'EXP', 'LOG'} within OPER }:
	y[n, 'EXP'] + y[(2*n+1), 'LOG'] <= 1;

subject to Red_LogExp {n in NLEAF: constr_redundancy >= 1 and 2*n+1 in NLEAF and {'EXP', 'LOG'} within OPER }:
	y[n, 'LOG'] + y[(2*n+1), 'EXP'] <= 1;

# [NEW] Aggregation of (2b, 3a, 3b) in [Cozad] no constaint unary operation, '- CST', or '/CST'
subject to RedNew_UDMCst {n in NLEAF: constr_redundancy >= 2 and constr_vdef_test_leaf == 0 }:
	y[2*n+1,'CST'] <= y[n,'PLUS'] + y[n,'MULT'];
	
# [NEW] No two minus on both node and its rchild, see A - (B - C) = A + (C - B)
subject to RedNew_DblMinus {n in NLEAF: constr_redundancy >= 3 }:
	y[n, 'PLUS'] + y[2*n+1, 'MINUS'] <= 1;

# [NEW] No two division on both node and its rchild, see A / (B / C) = A * (C / B)
subject to RedNew_DblDiv {n in NLEAF: constr_redundancy >= 3 }:
	y[n, 'MULT'] + y[2*n+1, 'DIV'] <= 1;

# [NEW] No (x,x) binary operation
#subject to RedNew_DblX {n in NLEAF, o in OPER_V: constr_redundancy >= 3}:
#	y[2*n, o] + y[(2*n+1), o] <= 1;

subject to Red_LeafModel {n in NLEAF: constr_redundancy >= 1 and constr_vdef_test_leaf == 1}:
	y[2*n+1, 'LEAF'] <= 1 - y[n,'PLUS'] - y[n,'MINUS'] ;


##### [NEW] Implication constraints avaiable when x in [L,U] where L<0 and U>0
param constr_implication default 0 ;
# (1-1) No log(x), sqrt(x), ()/x
subject to ImplicationLv1LogSqrtDiv {n in NLEAF: constr_implication >= 1}:
	sum {o in OPER_L} y[(2*n+1),o] + y[n,'LOG'] + y[n,'SQRT'] + y[n,'DIV'] <= 1;

# (1-2) No log(x), sqrt(x) when we know x != 0
subject to ImplicationLv1LogSqrt {n in NLEAF: constr_implication >= 1}:
	sum {o in OPER_L} y[(2*n+1),o] + y[n,'LOG'] + y[n,'SQRT'] <= 1;

# (2) No log(x^2), sqrt(x^2), ()/x^2 because x can be zero (we do not allow sqrt(0))  
subject to ImplicationLv2Zero {n in NODE: {2*n+1, 4*n+3} within NODE and constr_implication >= 2}:
	y[n,'LOG'] + y[n,'SQRT'] + y[n,'DIV']
	+ y[(2*n+1), 'SQR'] + sum {o in OPER_L} y[(2*n+1), o] 
	+ sum {o in OPER_L} y[4*n+3, o] <= 2; 

# (3) No log(C*x), sqrt(C*x) because x can be negative and positive
#     No ()/(C*x) because x can be zero 
subject to ImplicationLv2Sign1 {n in NODE: {2*n+1, 4*n+2, 4*n+3} within NODE and constr_implication >= 2}:
	y[n,'LOG'] + y[n,'SQRT']
	+ y[(2*n+1), 'MULT'] + sum {o in OPER_L} y[(2*n+1), o] 
	+ sum {o in OPER_V} y[4*n+2, o] 
	+ y[4*n+3, 'SQRT'] + y[4*n+3, 'EXP']
	+ sum {o in OPER_L} y[4*n+3, o]  <= 3;

subject to ImplicationLv2Sign2 {n in NODE: {2*n+1, 4*n+2, 4*n+3} within NODE and constr_implication >= 2}:
	y[n,'LOG'] + y[n,'SQRT']
	+ y[(2*n+1), 'MULT'] + sum {o in OPER_L} y[(2*n+1), o] 
	+ sum {o in OPER_V} y[4*n+3, o] 
	+ y[4*n+2, 'SQRT'] + y[4*n+2, 'EXP']
	+ sum {o in OPER_L} y[4*n+2, o]  <= 3;


##### CONSTR: Symmetry  #####
param constr_symmetry default 1 ;

# [Cozad] For the first data point, lch's value >= rch's value
subject to Symmetry {n in NLEAF: constr_symmetry >= 1}:
	v[first(OBS),2*n] - v[first(OBS),(2*n+1)] >= (v_lo - v_up) * (1 - y[n,'PLUS'] - y[n,'MULT']);




##### CONSTR: Value assignment constraints (defining v) #####
param constr_vdef default 1 ;
param constr_vdef_test default 0 ;

# When the node is None (inactive),
subject to ValueNoneUp {i in OBS, n in NODE: constr_vdef == 1}:
	v[i,n] <= v_up * sum{o in OPER} y[n,o];

subject to ValueNoneLo {i in OBS, n in NODE: constr_vdef == 1}:
	v[i,n] >= v_lo * sum{o in OPER} y[n,o];

# When the node is X,
subject to ValueXUp {i in OBS, n in NODE, o in OPER_V: constr_vdef == 1 and constr_vdef_test_leaf == 0}:
	x[i,o] - v[i,n] <= x_eps + (x[i,o] - v_lo) * (1 - y[n,o]);

subject to ValueXLo {i in OBS, n in NODE, o in OPER_V: constr_vdef == 1 and constr_vdef_test_leaf == 0}:
	x[i,o] - v[i,n] >= -x_eps + (x[i,o] - v_up) * (1 - y[n,o]);

# When the node is CST,
subject to ValueCstUpObs1 {n in NODE: constr_vdef == 1 and constr_vdef_test_leaf == 0}:
	v[first(OBS),n] <= cst_up * y[n,'CST'] + v_up * (1 - y[n,'CST']);

subject to ValueCstLoObs1 {n in NODE: constr_vdef == 1 and constr_vdef_test_leaf == 0}:
	v[first(OBS),n] >= cst_lo * y[n,'CST'] + v_lo * (1 - y[n,'CST']);

subject to ValueCstUp {i in OBS, n in NODE: ord(i) > 1 and constr_vdef >= 1 and constr_vdef_test_leaf == 0}:
	v[i,n] - v[first(OBS),n] <= (v_up - v_lo) * (1 - y[n,'CST']);

subject to ValueCstLo {i in OBS, n in NODE: ord(i) > 1 and constr_vdef >= 1 and constr_vdef_test_leaf == 0}:
	v[i,n] - v[first(OBS),n] >= (v_lo - v_up) * (1 - y[n,'CST']);


subject to ValueCstDomain {i in OBS, n in NODE: constr_vdef in {1,2,3} and constr_vdef_test_leaf == 0}:
	v[i,n]^2 >= cst_eps^2 * y[n,'CST'];


# LEAF test
subject to ValueLeafUp {i in OBS, n in NODE: constr_vdef >= 1 and constr_vdef_test_leaf == 1 }:
	v[i,n] - ( sum {o in OPER_V} (p[n,o] * x[i,o]) + c[n] ) <= v_up * (sum {(nn,oo) in YSET: oo in OPER_BU} y[nn,oo]) ;

subject to ValueLeafLo {i in OBS, n in NODE: constr_vdef >= 1 and constr_vdef_test_leaf == 1 }:
	v[i,n] - ( sum {o in OPER_V} (p[n,o] * x[i,o]) + c[n] ) >= v_lo * (sum {(nn,oo) in YSET: oo in OPER_BU} y[nn,oo]) ;

#subject to Leaf_DownCstX {n in NODE, o in OPER_V union {'CST'}: constr_vdef_test_leaf == 1}:
#	y[n,o] == 0 ;

#subject to DownLeaf {n in NODE: constr_vdef_test_leaf == 0}:
#	y[n,'LEAF'] == 0 ;

# [NEW] Alternative tight version of ValueNone, ValueCst, and ValueX1
subject to ValueLinearUp {i in OBS, n in NODE: constr_vdef >= 2}:
	v[i,n] <= ( if n in NLEAF then v_up * (sum{o in OPER_B} y[n,o]) else 0 )
		+ ( if n in NLEAF and 'SQR' in OPER then ( min(v_up, max(v_up^2, v_lo^2)) * y[n,'SQR'] ) else 0 )
		+ ( if n in NLEAF and 'SQRT' in OPER then ( min(v_up, sqrt(v_up)) * y[n,'SQRT'] ) else 0 )
		+ ( if n in NLEAF and 'EXP' in OPER then ( min(v_up, exp_scale * exp(v_up)) * y[n,'EXP'] ) else 0 )
		+ ( if n in NLEAF and 'LOG' in OPER then ( min(v_up, log_scale * log(v_up)) * y[n,'LOG'] ) else 0 )
		+ (if constr_vdef_test_leaf == 0 then cst_up * y[n,'CST'] else 0 )
		+ (if constr_vdef_test_leaf == 0 then sum {o in OPER_V} (x[i,o] * y[n,o]) else 0 )
		+ ( if 'LEAF' in OPER_L then v_up * y[n,'LEAF'] else 0 ) ;

subject to ValueLinearLo {i in OBS, n in NODE: constr_vdef >= 2}:
	v[i,n] >= ( if n in NLEAF then v_lo * (sum{o in OPER_B} y[n,o]) else 0 )
		+ ( if n in NLEAF and 'LOG' in OPER then ( max(v_lo, log_scale * log(log_eps)) * y[n,'LOG'] ) else 0 ) 
		+ ( if constr_vdef_test_leaf == 0 then cst_lo * y[n,'CST'] else 0 )
		+ ( if constr_vdef_test_leaf == 0 then sum {o in OPER_V} (x[i,o] * y[n,o]) else 0 )
		+ ( if 'LEAF' in OPER_L then v_lo * y[n,'LEAF'] else 0 );


subject to ValueQuadUp {i in OBS, n in NODE: constr_vdef == 2 and constr_vdef_test == 1}:
	v[i,n] <= ( if n in NLEAF then v_up * (y[n,'MULT'] + y[n,'DIV']) else 0 )
		+ ( if n in NLEAF then (v[i,2*n] + v[i,2*n+1]) * y[n,'PLUS'] else 0 )
		+ ( if n in NLEAF then (v[i,2*n] - v[i,2*n+1]) * y[n,'MINUS'] else 0 )
		+ ( if n in NLEAF and 'SQR' in OPER then ( min(v_up, max(v_up^2, v_lo^2)) * y[n,'SQR'] ) else 0 )
		+ ( if n in NLEAF and 'SQRT' in OPER then ( min(v_up, sqrt(v_up)) * y[n,'SQRT'] ) else 0 )
		+ ( if n in NLEAF and 'EXP' in OPER then ( min(v_up, exp_scale * exp(v_up)) * y[n,'EXP'] ) else 0 )
		+ ( if n in NLEAF and 'LOG' in OPER then ( min(v_up, log_scale * log(v_up)) * y[n,'LOG'] ) else 0 )
		+ ( if ord(i) == 1 then cst_up * y[n,'CST'] else v[first(OBS),n] * y[n,'CST'] )
		+ sum {o in OPER_V} (x[i,o] * y[n,o]);

subject to ValueQuadLo {i in OBS, n in NODE: constr_vdef == 2 and constr_vdef_test == 1}:
	v[i,n] >= ( if n in NLEAF then v_lo * (y[n,'MULT'] + y[n,'DIV']) else 0 )
		+ ( if n in NLEAF then (v[i,2*n] + v[i,2*n+1]) * y[n,'PLUS'] else 0 )
		+ ( if n in NLEAF then (v[i,2*n] - v[i,2*n+1]) * y[n,'MINUS'] else 0 )
		+ ( if n in NLEAF and 'LOG' in OPER then ( max(v_lo, log_scale * log(log_eps)) * y[n,'LOG'] ) else 0 ) 
		+ ( if ord(i) == 1 then cst_lo * y[n,'CST'] else v[first(OBS),n] * y[n,'CST'] )
		+ sum {o in OPER_V} (x[i,o] * y[n,o]) ;

subject to ValueBoundRchUp {i in OBS, n in NLEAF: constr_vdef >= 2 and ('SQR' in OPER or 'EXP' in OPER) }:
	v[i,2*n+1] <= v_up * (sum {o in OPER_BU diff {'SQR', 'EXP'}} y[n,o])
			+ ( if 'SQR' in OPER then min(v_up, sqrt(v_up)) * y[n,'SQR'] else 0 )
			+ ( if 'EXP' in OPER then min(v_up, log(v_up))  * y[n,'EXP'] else 0 ) ;

subject to ValueBoundRchLo {i in OBS, n in NLEAF: constr_vdef >= 2 and ('SQRT' in OPER or 'LOG' in OPER) }:
	v[i,2*n+1] >= v_lo * (sum {o in OPER_BU diff {'SQRT', 'LOG'}} y[n,o])
			+ ( if 'SQRT' in OPER then sqrt_eps * y[n,'SQRT'] else 0 )
			+ ( if 'LOG' in OPER then log_eps * y[n,'LOG'] else 0 );

# When the node is PLUS,
subject to ValuePlusUp {i in OBS, n in NLEAF: constr_vdef in {1,2} and constr_vdef_test == 0}:
	(v[i,2*n] + v[i,(2*n+1)]) - v[i,n] <= (2*v_up - v_lo) * (1 - y[n,'PLUS']);

subject to ValuePlusLo {i in OBS, n in NLEAF: constr_vdef in {1,2} and constr_vdef_test == 0}:
	(v[i,2*n] + v[i,(2*n+1)]) - v[i,n] >= (2*v_lo - v_up) * (1 - y[n,'PLUS']);

# [NEW] tight version of When the node is PLUS,
subject to ValuePlusUp2 {i in OBS, n in NLEAF: constr_vdef == 3}:
	v[i,n] - (v[i,2*n] + v[i,(2*n+1)]) <= (v_up - 2*v_lo) * (sum {o in OPER_B diff {'PLUS'}} y[n,o])
						+ (v_up - v_lo) * (sum {o in OPER_U} y[n,o])
						+ cst_up * y[n,'CST']
						+ sum {o in OPER_V} (x[i,o] * y[n,o]);

subject to ValuePlusLo2 {i in OBS, n in NLEAF: constr_vdef == 3}:
	v[i,n] - (v[i,2*n] + v[i,(2*n+1)]) >= (v_lo - 2*v_up) * (sum {o in OPER_B diff {'PLUS'}} y[n,o])
						+ (v_lo - v_up) * (sum {o in OPER_U} y[n,o])
						+ cst_lo * y[n,'CST']
						+ sum {o in OPER_V} (x[i,o] * y[n,o]);


# When the node is MINUS,
subject to ValueMinusUp {i in OBS, n in NLEAF: constr_vdef in {1,2} and constr_vdef_test == 0}:
	(v[i,2*n] - v[i,(2*n+1)]) - v[i,n] <= (v_up - 2*v_lo) * (1 - y[n,'MINUS']);

subject to ValueMinusLo {i in OBS, n in NLEAF: constr_vdef in {1,2} and constr_vdef_test == 0}:
	(v[i,2*n] - v[i,(2*n+1)]) - v[i,n] >= (v_lo - 2*v_up) * (1 - y[n,'MINUS']);

# [NEW] tight version of When the node is MINUS,
subject to ValueMinusUp2 {i in OBS, n in NLEAF: constr_vdef == 3}:
	v[i,n] - (v[i,2*n] - v[i,(2*n+1)]) <= (2*v_up - v_lo) * (sum {o in OPER_B diff {'MINUS'}} y[n,o])
						+ 2*v_up * (sum {o in OPER_U} y[n,o])
						+ cst_up * y[n,'CST']
						+ sum {o in OPER_V} (x[i,o] * y[n,o]);

subject to ValueMinusLo2 {i in OBS, n in NLEAF: constr_vdef == 3}:
	v[i,n] - (v[i,2*n] - v[i,(2*n+1)]) >= (2*v_lo - v_up) * (sum {o in OPER_B diff {'MINUS'}} y[n,o])
						+ 2*v_lo * (sum {o in OPER_U} y[n,o])
						+ cst_lo * y[n,'CST']
						+ sum {o in OPER_V} (x[i,o] * y[n,o]);


# When the node is MULT,
subject to ValueMultUp {i in OBS, n in NLEAF: constr_vdef in {1,2}}:
	(v[i,2*n] * v[i,(2*n+1)]) - v[i,n] <= (max(v_up^2, v_lo^2, v_up*v_lo) - v_lo) * (1 - y[n,'MULT']);

subject to ValueMultLo {i in OBS, n in NLEAF: constr_vdef in {1,2}}:
	(v[i,2*n] * v[i,(2*n+1)]) - v[i,n] >= (min(v_up^2, v_lo^2, v_up*v_lo) - v_up) * (1 - y[n,'MULT']);

# [NEW] tight version of When the nodeis MULT,
subject to ValueMultUp2 {i in OBS, n in NLEAF: constr_vdef == 3}:
	v[i,n] - (v[i,2*n] * v[i,(2*n+1)]) <= (v_up - min(v_up^2, v_lo^2, v_up*v_lo)) * (sum {o in OPER_B diff {'MULT'}} y[n,o])
					      	+ v_up * (sum {o in OPER_U} y[n,o])
						+ cst_up * y[n,'CST']
						+ sum {o in OPER_V} (x[i,o] * y[n,o]) ;

subject to ValueMultLo2 {i in OBS, n in NLEAF: constr_vdef == 3}:
	v[i,n] - (v[i,2*n] * v[i,(2*n+1)]) >= (v_lo - max(v_up^2, v_lo^2, v_up*v_lo)) * (sum {o in OPER_B diff {'MULT'}} y[n,o])
						+ v_lo * (sum {o in OPER_U} y[n,o])
						+ cst_lo * y[n,'CST']
						+ sum {o in OPER_V} (x[i,o] * y[n,o]);


subject to ValueMultDomainLch {i in OBS, n in NLEAF: constr_vdef in {1,2,3}}:
	v[i,(2*n)]^2 >= mult_eps^2 * y[n,'MULT'];

subject to ValueMultDomainRch {i in OBS, n in NLEAF: constr_vdef in {1,2,3}}:
	v[i,(2*n+1)]^2 >= mult_eps^2 * y[n,'MULT'];

# When the node is DIV,
subject to ValueDivUp {i in OBS, n in NLEAF: constr_vdef in {1,2}}:
	v[i,2*n] - v[i,(2*n+1)] * v[i,n] <= (v_up - min(v_up^2, v_lo^2, v_up*v_lo)) * (1 - y[n,'DIV']);

subject to ValueDivLo {i in OBS, n in NLEAF: constr_vdef in {1,2}}:
	v[i,2*n] - v[i,(2*n+1)] * v[i,n] >= (v_lo - max(v_up^2, v_lo^2, v_up*v_lo)) * (1 - y[n,'DIV']);

subject to ValueDivDomain {i in OBS, n in NLEAF: constr_vdef in {1,2,3}}:
	v[i,(2*n+1)]^2 >= div_eps^2 * y[n,'DIV'];

# [NEW] tight version of When the node is DIV (when OPER_L, LHS becomes zero)
subject to ValueDivUp2 {i in OBS, n in NLEAF: constr_vdef in {3}}:
	v[i,2*n] - v[i,(2*n+1)] * v[i,n] <= (v_up - min(v_up^2, v_lo^2, v_up*v_lo)) 	* (sum {o in OPER_B diff {'DIV'}} y[n,o])
						+ max(-v_up^3, -v_lo^3) 		* y[n,'SQR']
						- v_lo * exp_scale * exp(v_lo) 		* y[n,'EXP']
						- (log_eps * v_lo) 			* y[n,'LOG'] ;

subject to ValueDivLo2 {i in OBS, n in NLEAF: constr_vdef in {3}}:
	v[i,2*n] - v[i,(2*n+1)] * v[i,n] >= (v_lo - max(v_up^2, v_lo^2, v_up*v_lo)) 	* (sum {o in OPER_B diff {'DIV'}} y[n,o])
						+ min(-v_up^3, -v_lo^3) 		* y[n,'SQR']
						+ - max(0,v_up)^(3/2) 			* y[n,'SQRT']
						+ max(- v_up^2, - v_up * exp_scale * exp(v_up))	* y[n,'EXP']
					        + max(- v_up^2, - v_up * log_scale * log(v_up)) * y[n,'LOG'] ;

	
# When the node is SQR,
subject to ValueSqrUp {i in OBS, n in NLEAF: constr_vdef in {1,2} and 'SQR' in OPER }:
	v[i,(2*n+1)]^2 - v[i,n] <= (max(v_up^2, v_lo^2) - v_lo) * (1 - y[n,'SQR']);

subject to ValueSqrLo {i in OBS, n in NLEAF: constr_vdef in {1,2} and 'SQR' in OPER }:
	v[i,(2*n+1)]^2 - v[i,n] >= -v_up  * (1 - y[n,'SQR']);

# [NEW] tight version of SQR
subject to ValueSqrUp2 {i in OBS, n in NLEAF: constr_vdef in {3} and 'SQR' in OPER }:
	v[i,n] - v[i,(2*n+1)]^2 <= v_up * (sum {o in OPER_BU diff {'SQR'}} y[n,o])
					+ cst_up * y[n,'CST']
					+ sum {o in OPER_V} (x[i,o] * y[n,o]) ;

subject to ValueSqrLo2 {i in OBS, n in NLEAF: constr_vdef in {3} and 'SQR' in OPER }:
	v[i,n] - v[i,(2*n+1)]^2 >= min(-v_up^2, -v_lo^2) * (sum {o in OPER_BU diff {'SQR'}} y[n,o])
					+ cst_lo * y[n,'CST']
					+ sum {o in OPER_V} (x[i,o] * y[n,o]) ;			   


# When the node is SQRT,
subject to ValueSqrtUp {i in OBS, n in NLEAF: 'SQRT' in OPER and constr_vdef in {1,2} }:
	v[i,(2*n+1)] - v[i,n]^2 <= v_up * (1 - y[n,'SQRT']);

subject to ValueSqrtLo {i in OBS, n in NLEAF: 'SQRT' in OPER and constr_vdef in {1,2} }:
	v[i,(2*n+1)] - v[i,n]^2 >= (v_lo - max(v_up^2, v_lo^2)) * (1 - y[n,'SQRT']);

# [NEW] tight version of SQRT
subject to ValueSqrtUp2 {i in OBS, n in NLEAF: 'SQRT' in OPER and constr_vdef in {3} }:
	v[i,n]^2 - v[i,(2*n+1)] <= (max(v_up^2, v_lo^2) - v_lo) * (sum {o in OPER_BU diff {'SQRT'}} y[n,o])
 					+ cst_up^2 * y[n,'CST']
					+ sum {o in OPER_V} (x[i,o]^2 * y[n,o]);

subject to ValueSqrtLo2 {i in OBS, n in NLEAF: 'SQRT' in OPER and constr_vdef in {3} }:
	v[i,n]^2 - v[i,(2*n+1)] >= (min(v_up^2, v_lo^2, v_up*v_lo) - v_up) * (sum {o in OPER_BU diff {'SQRT'}} y[n,o])
					+ 0 * y[n,'CST']
					+ sum {o in OPER_V} (x[i,o]^2 * y[n,o]);

# Domain restriction on SQRT
subject to ValueSqrtDomain {i in OBS, n in NLEAF: 'SQRT' in OPER and constr_vdef in {1,2,3}}:
	v[i,(2*n+1)] >= sqrt_eps * y[n,'SQRT'] + v_lo * (1 - y[n,'SQRT']);


# When the node is EXP,
subject to ValueExpUp {i in OBS, n in NLEAF: 'EXP' in OPER and constr_vdef in {1,2} }:
	v[i,n] - exp_scale * exp(v[i,(2*n+1)])
		<= (exp_scale * exp(v_up) - v_lo) * (1 - y[n,'EXP']);

subject to ValueExpLo {i in OBS, n in NLEAF: 'EXP' in OPER and constr_vdef in {1,2} }:
	v[i,n] - exp_scale * exp(v[i,(2*n+1)])
		>= (exp_scale * exp(v_lo) - v_up) * (1 - y[n,'EXP']);

# [NEW] Tight version of EXP
subject to ValueExpUp2 {i in OBS, n in NLEAF: 'EXP' in OPER and constr_vdef in {3} }:
	v[i,n] - exp_scale * exp(v[i,(2*n+1)])
		<= -1 		* (1 - sum {o in OPER} y[n,o])	# None
		+ (v_up - exp_scale * exp(v_lo)) * (sum {o in OPER_BU diff {'EXP'}} y[n,o])
		+ (v_up - 1) 	* y[n,'CST']
		+ sum {o in OPER_V} ( (x[i,o] - 1) * y[n,o] ) ;

subject to ValueExpLo2 {i in OBS, n in NLEAF: 'EXP' in OPER and constr_vdef in {3} }:
	v[i,n] - exp_scale * exp(v[i,(2*n+1)]) 
		>= -1					* (1 - sum {o in OPER} y[n,o])	# None
		+ (v_lo - exp_scale * exp(v_up)) 	* (sum {o in OPER_BU diff {'EXP'}} y[n,o])
		+ (v_lo - 1) 				* y[n,'CST']
		+ sum {o in OPER_V} ( (x[i,o] - 1) * y[n,o] );

# When the node is LOG,
subject to ValueLogUp {i in OBS, n in NLEAF: 'LOG' in OPER and constr_vdef in {1,2} }:
	exp(v[i,n] / log_scale) - v[i,(2*n+1)]  
		<= (exp(v_up / log_scale) - v_lo) * (1 - y[n,'LOG']);

subject to ValueLogLo {i in OBS, n in NLEAF: 'LOG' in OPER and constr_vdef in {1,2} }:
	exp(v[i,n] / log_scale) - v[i,(2*n+1)]  
		>= (exp(v_lo / log_scale) - v_up) * (1 - y[n,'LOG']);

# [NEW] tight version of LOG
subject to ValueLogUp2 {i in OBS, n in NLEAF: 'LOG' in OPER and constr_vdef in {3} }:
	exp(v[i,n] / log_scale) - v[i,(2*n+1)]  
		<= 1					* (1 - sum {o in OPER} y[n,o])	# None
		+ exp(cst_up / log_scale)		* y[n,'CST']
		+ sum {o in OPER_V} ( exp(x[i,o] / log_scale) * y[n,o] )
		+ (exp(v_up / log_scale) - v_lo)	* (sum {o in OPER_BU diff {'LOG'}} y[n,o]);
	
subject to ValueLogLo2 {i in OBS, n in NLEAF: 'LOG' in OPER and constr_vdef in {3} }:
	exp(v[i,n] / log_scale) - v[i,(2*n+1)]  
		>= 1					* (1 - sum {o in OPER} y[n,o])	# None
		+ exp(cst_lo / log_scale)		* y[n,'CST']
		+ sum {o in OPER_V} ( exp(x[i,o] / log_scale) * y[n,o] )
		+ (exp(v_lo / log_scale) - v_up)	* (sum{o in OPER_BU diff {'LOG'}} y[n,o]);

# Domain restriction on LOG
subject to ValueLogDomain {i in OBS, n in NLEAF: 'LOG' in OPER and constr_vdef in {1,2,3} }:
	v[i,(2*n+1)] >= log_eps * y[n,'LOG'] + v_lo * (1 - y[n,'LOG']);





##### Alternative node value constraints #####
# CST X1 PLUS MINUS MULT DIV SQR SQRT EXP LOG
#var v_div 	{i in OBS, n in NLEAF};
#var v_sqrt 	{i in OBS, n in NLEAF};
#var v_log 	{i in OBS, n in NLEAF};

#param formulation_NV default 0 ;

#subject to NodeValueNleaf {i in OBS, n in NLEAF: formulation_NV}:
#	v[i,n] == y[n,'CST'] * (v[1,n])
#		+ sum {o in OPER_V} (y[n,o] * x[i,o])
#		+ y[n,'PLUS'] * (v[i,2*n] + v[i,(2*n+1)])
#		+ y[n,'MINUS'] * (v[i,2*n] - v[i,(2*n+1)])
#		+ y[n,'MULT'] * (v[i,2*n] * v[i,(2*n+1)])
#		+ y[n,'DIV'] * v_div[i,n]
#		+ y[n,'SQR'] * (v[i,(2*n+1)])^2
#		+ y[n,'SQRT'] * v_sqrt[i,n]
#		+ y[n,'EXP'] * (exp(log(v_up) / v_up * v[i,(2*n+1)]))
#		+ y[n,'LOG'] * v_log[i,n] ;

#subject to NodeValueLeaf {i in OBS, n in LEAF: formulation_NV}:
#	v[i,n] == y[n,'CST'] * (v[1,n]) + sum {o in OPER_V} (y[n,o] * x[i,o]);

#subject to NodeValueDiv {i in OBS, n in NLEAF: formulation_NV}:
#	v_div[i,n] * v[i,(2*n+1)] == v[i,2*n];

#subject to NodeValueSqrt {i in OBS, n in NLEAF: formulation_NV}:
#	v_sqrt[i,n]^2 == v[i,(2*n+1)];

#subject to NodeValueLog {i in OBS, n in NLEAF: formulation_NV}:
#	exp(log(v_up) / v_up * v_log[i,n]) == v[i,(2*n+1)];










