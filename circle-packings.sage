'''
Basic data structures for circle packings:
==========================================

The code computes a circle packing numerically from a given planar graph.
It is used to search for algebraic relations in circle packings.

Author: Benjamin Matschke

License: Creative Commons 4.0 BY-NC.
'''

from mpmath import pslq #in particular pslq

def randomDict(keys):
	d = {};
	for k in keys:
		d[k] = random()+0.5;
	return d;

def monomials(xList, degree):
	result = [];
	n = len(xList);
	for A in unordered_tuples(range(n),degree):
		monomial = prod([xList[i] for i in A]);
		result.append(monomial);
	return result;

class CirclePacking:

	def __init__(self,G,radiiOfFixedNodes=None,center=vector([0,0]),angle=0,verbose=False):
		'''
		G is the graph of a triangulation of the sphere, with a prescribed midpoints. 
		radiiOfBoundaryNodes is a dictionary, whose keys are the boundary nodes,
		   and whose values are their prescribed radii.
		'''
		
		self.G = G;
		self.center = center; #for midpoints in R^2
		self.angle = angle; #for midpoints in R^2
		if not G.is_planar(set_embedding = True):
			print "ValueError: Given graph G is non-planar.";
		self.neighborsOfVertex = G.get_embedding();

		if radiiOfFixedNodes == None:
			v0 = G.vertices()[0];
			v1 = self.neighborsOfVertex[v0][0];
			v2 = self.neighborsOfVertex[v0][-1];
			radiiOfFixedNodes = randomDict([v0,v1,v2]);
		self.radiiOfFixedNodes = radiiOfFixedNodes;

		self.RR = RealField(200);
		eps = self.RR(10)^(-30);
		delta = self.RR(10)^(-50);

		self.fixedNodes = radiiOfFixedNodes.keys();
		self.freeNodes = [u for u in self.G.vertices() if u not in self.fixedNodes];

		self.R = {};
		for v,r in radiiOfFixedNodes.iteritems():
			self.R[v] = self.RR(r);
		maxR = max(radiiOfFixedNodes.values());
		for u in self.freeNodes:
			self.R[u] = maxR;

		loop = 0;
		c = eps + 1;
		lamda = self.RR(-1);
		flag = False;
		while c>eps:
			loop += 1;
			if verbose:
				print "Loop",loop,":";

			#Step (a):
			c0 = c;
			lamda0 = lamda;
			flag0 = flag;
			R0 = copy(self.R);
			#print "R0:",R0;

			#Step (b):
			newR = self.newRadii();
			self.R = newR;
			defect = self.totalAngleDefects();
			c += defect;
			if verbose:
				print "defect:",defect;
				print "c =",c;
			
			#Step (c):
			c = sqrt(c); #This c is roughly the l2-norm of the angle defect vector (and it remembers some history of c).
			lamda = c/c0;
			flag = True;
			#flag = False; #debug: No acceleration
			if verbose:
				print "lamda =",lamda;

			#Step (d):
			if flag0 and lamda<1:
				c = lamda*c;
				if abs(lamda-lamda0)<delta:
					lamda = lamda/(1-lamda);
					if verbose:
						print "Do SUPER acceleration!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1!"
				else:
					if verbose:
						print "Do simple acceleration."
				lamdaStar = max([-self.R[v]/(self.R[v]-R0[v]) for v in self.R.keys() if self.R[v]<R0[v]]+[0]);
				lamda = min(lamda,self.RR(0.005)*lamdaStar);
				for v in self.R.keys():
					self.R[v] = max(eps,self.R[v] + lamda*(self.R[v] - R0[v]));
				flag = False;

			if defect < eps:
				break;
			if loop >= 500:
				break;

		self.curvatures = dict([(u,1/self.R[u]) for u in G.vertices()]);

		self.embedIntoThePlane();

		n = G.num_verts();
		self.gramMatrix = matrix(RR,n,n);
		for i in range(n):
			u = G.vertices()[i];
			C1 = self.circles[u];
			for j in range(n):
				v = G.vertices()[j];
				C2 = self.circles[v];
				self.gramMatrix[i,j] = C1.pairingWith(C2);

		if verbose:
			print "R =",self.R;
			print "Free nodes:",self.freeNodes;
			print "NeighborsOfVertex:",self.neighborsOfVertex;

			print "Done.";
			print "Radii:",self.R;
			print "Radii of free circles as continued fractions:"
			for v in self.freeNodes:
				print v,":",continued_fraction_list(self.R[v]);

			print "Curvatures:",self.curvatures;
			print "Curvatures of free circles as continued fractions:"
			for v in self.freeNodes:
				print v,":",continued_fraction_list(self.curvatures[v]);
			
	def angleAtX(self,x,y,z):
		'''
		Given are three pairwise tangent circles with radii x,y,z.
		Their midpoints form a triangle.
		The output is the interior angle of that triangle
		  at the vertex corresponding to x.
		'''
		
		alpha = 2*arcsin(sqrt(y/(x+y)*z/(x+z)));
		alpha = self.RR(alpha);
	
		#Should yield the same as:
		#alpha = arccos(((x+y)^2+(x+z)^2-(y+z)^2)/(2*(x+y)*(x+z)));
		
		return alpha;

	def angleSumAtU(self,u):
		theta = 0;
		neighborsOfU = self.neighborsOfVertex[u];
		R = self.R;
		#print "self.R =",self.R;
		#print "R:",R;
		for i in range(len(neighborsOfU)):
			v = neighborsOfU[i];
			w = neighborsOfU[(i+1) % len(neighborsOfU)];
			#print "u,w:",u,w;
			#print "R[u] =",R[u];
			alpha = self.angleAtX(R[u],R[v],R[w]);
			theta += alpha;
		return theta;
	
	def newRadiusAtU(self,u,verbose=False):
		#print "u =",u;
		A = 2*pi; #The desired angle sum at u.
		neighborsOfU = self.neighborsOfVertex[u];
		k = len(neighborsOfU);
		theta = self.angleSumAtU(u);
		beta = sin(theta/(2*k));
		delta = sin(A/(2*k));
		radius_hat = beta/(1-beta) * self.R[u];
		radius_new = (1-delta)/delta * radius_hat;
		if verbose:
			print "k =",k;
			print "theta =",theta;
			print "beta =",beta;
			print "delta =",delta;
			print "radius_hat =",radius_hat;
			print "radius_new =",radius_new;
		return self.RR(radius_new);

	def newRadii(self):
		result = copy(self.R);
		for u in self.freeNodes:
			result[u] = self.newRadiusAtU(u);
		return result;

	def totalAngleDefects(self,exponent=2):
		A = 2*pi;
		result = 0;
		for u in self.freeNodes:
			defectAtU = abs(A - self.angleSumAtU(u));
			result += self.RR(defectAtU)^exponent;
		return result;

	def embedIntoThePlane(self,verbose=False):
		G = self.G;
		R = self.R;
		neighborsOfVertex = self.neighborsOfVertex;
		midpoints = {};
		v0 = G.vertices()[0];
		v1 = G.edges_incident(v0)[0][1];
		midpoints[v0] = self.center;
		midpoints[v1] = self.center + (R[v0]+R[v1]) * vector([cos(self.angle),sin(self.angle)]);

		verticesToVisitNext = [v1,v0];
		while verticesToVisitNext != []:
			v = verticesToVisitNext.pop();
			#print "visiting v=",v;
			neighbors = neighborsOfVertex[v];
			#print "neighbors",neighbors;
			for i in range(len(neighbors)):
				if neighbors[i] in midpoints.keys():
					i0 = i;
			for i in range(i0+1,len(neighbors)) + range(i0):
				w0 = neighbors[i-1]; #its position in the plane is known
				w = neighbors[i]; #its position in the plane is to be determined

				if not midpoints.has_key(w):
					angle = self.angleAtX(R[v],R[w0],R[w]);
					factor = self.RR((R[v]+R[w])/(R[v]+R[w0]));
					midpoints[w] = self.rotateAndScale2dVector(midpoints[w0],angle,midpoints[v],factor);
					verticesToVisitNext.append(w);
					if verbose:
						print "obtain",w,"from rotating",w0,"around",v,"by angle",RR(angle*180/pi);
		self.midpoints = midpoints;
		self.circles = dict([(v,CircleInPlane(w=None,M=PointOnPlane(midpoints[v]),radius=R[v])) for v in G.vertices()]);
		return midpoints;
					
	def rotateAndScale2dVector(self,v,angle,center=vector([0,0]),factor=1):
		v0 = v-center;
		M = matrix([[cos(angle),-sin(angle)],[sin(angle),cos(angle)]]);
		result = center + factor * M * v0;
		#print "v0:",v0;
		#print "M:",M;
		#print "result:",result;
		#print "center:",center;
		return result;

	def show(self,verbose=False):
		g = Graphics();
		#midpoints = self.embedIntoThePlane(); #was done already during __init__()
		midpoints = self.midpoints;
		if verbose:
			print "midpoints:",midpoints;
		for (v,m) in midpoints.iteritems():
			#print "v,m:",v,m;
			radius = self.R[v];
			c = circle(m,radius);
			g += c;
			g += text(str(v),m,fontsize = max(10,min(4,radius*g.SHOW_OPTIONS['dpi'])));
		g.show(axes=False);

	def relationsAmongCurvatures(self):
		cs = self.curvatures.values();
		print "Curvatures:",cs;
		print "Relations among curvatures:";
		print pslq(cs,tol=0.00000001);

def findLinearRelationsAmongCurvatures(G):
	n = G.num_verts();
	M = matrix(ZZ,2*n,n);
	precision = 10;
	factor = 10^precision;
	for i in range(n):
		M[i,i] = 1;
	for i in range(n):
		cp = CirclePacking(G);
		curvatures = [cp.curvatures[u] for u in G.vertices()];
		for j in range(n):
			M[n+i,j] = round(factor*curvatures[j]);
	#print "cols of M =";
	#for col in M.columns():
	#	print col;
	Mred = M.transpose().LLL().transpose(); #LLL reduces columns of M
	print "LLL reduced columns:";
	for col in Mred.columns():
		print col;

def findQuadraticRelationsAmongCurvatures(G):
	n = G.num_verts();
	m = len(monomials(range(n),1) + monomials(range(n),2));
	M = matrix(ZZ,2*m,m);
	precision = 10;
	factor = 10^precision;
	for i in range(m):
		M[i,i] = 1;
	for i in range(m):
		cp = CirclePacking(G);
		curvatures = [cp.curvatures[u] for u in G.vertices()];
		row = monomials(curvatures,1) + monomials(curvatures,2);
		for j in range(m):
			M[m+i,j] = round(factor*row[j]);
	print "M =";
	for col in M.columns():
		print col;
	Mred = M.transpose().LLL().transpose(); #LLL reduces columns of M
	print "LLL reduced columns:";
	for col in Mred.columns():
		print col;

	varNames = [];
	for i in range(n):
		varNames.append('x'+str(i));
	R = PolynomialRing(QQ,varNames);
	monoms = monomials(R.gens(),1) + monomials(R.gens(),2);
	relations = [];
	for col in Mred.columns():
		if vector(col).norm() <= 10 * m:
			relation = sum([col[i]*monoms[i] for i in range(m)]);
			relations.append(relation);
			print "Relation: 0 =",relation;

	I = Ideal(relations);
	print "Ideal of relations:",I;
	GB = I.groebner_basis();
	print "Groebner basis of I:",GB;
	return I;

class CirclePackingVariety:

	def __init__(self,G):
		self.G = G;
		self.variableNames = [];
		self.vertexToIndex = {};
		i = 0;
		for v in G.vertices():
			rName = "r"+str(v);
			xName = "x"+str(v);
			yName = "y"+str(v);
			if i in [0,1]:
				self.variableNames += [rName];
			else:
				self.variableNames += [rName,xName,yName];
			self.vertexToIndex[v] = i;
			i += 1;
		self.R = PolynomialRing(QQ,self.variableNames);
		print "R =",self.R;
		self.n = G.num_verts();
		self.rs = [self.R.gen(0),self.R.gen(1)];
		self.xs = [self.R(0),self.R(1)];
		self.ys = [self.R(0),self.R(0)];
		for i in range(2,self.n):
			self.rs.append(self.R.gen(3*i-4));
			self.xs.append(self.R.gen(3*i-3));
			self.ys.append(self.R.gen(3*i-2));
		print "rs:",self.rs;
		print "xs:",self.xs;
		print "ys:",self.ys;
		self.relations = [];
		for e in G.edges():
			i = self.vertexToIndex[e[0]];
			j = self.vertexToIndex[e[1]];
			eRelation = (self.xs[i]-self.xs[j])^2 + (self.ys[i]-self.ys[j])^2 - (self.rs[i]+self.rs[j])^2;
			self.relations.append(eRelation);
			print -eRelation;
		self.A = Ideal(self.relations);
		print "A = ",self.A;

	def checkLinearRelation(self,v1,v2, w1,w2):
		i1 = self.vertexToIndex[v1];
		i2 = self.vertexToIndex[v2];
		j1 = self.vertexToIndex[w1];
		j2 = self.vertexToIndex[w2];
		r1 = self.rs[i1];
		r2 = self.rs[i2];
		r3 = self.rs[j1];
		r4 = self.rs[j2];
		#p = self.R(r1*r2*r3*r4*(1/r1+1/r2-1/r3-1/r4));
		p = r2*r3*r4 + r1*r3*r4 - r1*r2*r4 - r1*r2*r3;
		relationHolds = p in self.A;
		print "relation holds:",relationHolds;
		return relationHolds;

def findRelationsInCirclePacking(G,degrees=[1],whichValues=lambda cp: [cp.curvatures[u] for u in cp.G.vertices()],nSamples=10,varNames=None):
	precision = 10;
	factor = 10^precision;
	rows = [];
	print "Run through samples:",;
	for i in range(nSamples):
		print i,;
		sys.stdout.flush();
		cp = CirclePacking(G,center=vector([random(),random()]),angle=random());
		values = whichValues(cp);
		n = len(values); #should be the same for all samples
		row = [];
		for deg in degrees:
			row += monomials(values,deg);
		rowApprox = [round(factor*mon) for mon in row];
		rows.append(rowApprox);
	m = len(rows[0]); #number of monomials per sample
	rows = identity_matrix(ZZ,m).rows() + rows;
	M = matrix(ZZ,rows);
	#print "M =";
	#for col in M.columns():
	#	print col;
	Mred = M.transpose().LLL().transpose(); #LLL reduces columns of M
	print "LLL reduced columns:";
	for col in Mred.columns():
		print col;
	sys.stdout.flush();
	
	return;

	if varNames == None:
		varNames = [];
		for i in range(n):
			varNames.append('x'+str(i));
	R = PolynomialRing(QQ,varNames);
	monoms = [];
	for deg in degrees:
		monoms += monomials(R.gens(),deg);
	print "monoms:",monoms;
	sys.stdout.flush();
	relations = [];
	for col in Mred.columns():
		if vector(col).norm() <= 10 * m:
			relation = sum([col[i]*monoms[i] for i in range(m)]);
			relations.append(relation);
			if relation.is_prime():
				print "Relation: 0 =",relation;
			sys.stdout.flush();

	I = Ideal(relations);
	print "Ideal of relations:",I;
	sys.stdout.flush();
	GB = I.groebner_basis();
	print "Groebner basis of I:",GB;
	sys.stdout.flush();
	return I;

class PointOnPlane:
	def __init__(self,V):
		self.V = V;
		self.X = V[0];
		self.Y = V[1];

	def __str__(self):
		return "("+str(self.X)+","+str(self.Y)+")";
	
	def toSphere(self):
		X = self.X;
		Y = self.Y;
		denom = 1 + X^2 + Y^2;
		x = 2*X / denom;
		y = 2*Y / denom;
		z = -(denom-2)/denom; #steriographic projection from SOUTH pole (0,0,-1)!
		return PointOnSphere(vector([x,y,z]));

	def negative(self):
		return self.toSphere().negative().toPlane();

	def inversion(self):
		rSq = self.X^2 + self.Y^2;
		return PointOnPlane((1/rSq) * self.V);

	def translation(self,W):
		return PointOnPlane(self.V + W);

	def moebiusTransformation(self,W1,W2):
		return self.inversion().translation(W1).inversion().translation(W2);

	def normSquared(self):
		return self.X^2 + self.Y^2;

class PointOnSphere:
	def __init__(self,v):
		self.v = v;
		self.x = v[0];
		self.y = v[1];
		self.z = v[2];

	def __str__(self):
		return "("+str(self.x)+","+str(self.y)+str(self.z)+")";

	def toPlane(self):
		X = simplify(self.x/(1+self.z));
		Y = simplify(self.y/(1+self.z));
		return PointOnPlane(vector([X,Y]));

	def negative(self):
		return PointOnSphere(-self.v);

	def moebiusTransformation(self,W1,W2):
		return self.toPlane().moebiusTransformation(W1,W2).toSphere();

class CircleInPlane:

	def __init__(self,w=None,M=None,radius=None):
		'''
		Either curvature coordinates w must be given, or
		a PointOnPlane M and a positive real radius.
		'''
		
		if M!=None and radius!=None:
			self.M = M; #A PointOnPlane
			self.radius = radius;
			if w==None:
				c = 1/radius;
				cBar = c*M.normSquared() - radius;
				self.w = vector([cBar,c,c*M.X,c*M.Y]);
			else:
				self.w = w;
		elif w!=None:
			self.w = w;
			if radius==None:
				self.radius = 1/w[1];
			else:
				self.radius = radius;
			if M==None:
				c = w[1];
				self.M = PointOnPlane(vector([w[2]/c,w[3]/c]));
			else:
				self.M = M;
		else:
			print "Error: Not enough parameters to determine the given circle!";
		

	def __str__(self):
		return "Circle: M = "+str(self.M)+", r = "+str(self.radius)+".";

	def toSphere(self):
		w = self.w;
		wPlus = vector([(w[0]+w[1])/2,(w[0]-w[1])/2,w[2],w[3]]);
		return CircleInSphere(wPlus=wPlus);

	def negative(self):
		return self.toSphere().negative().toPlane();

	def translation(self,v):
		w = self.w;
		c = w[1];
		new_w = vector([w[0],w[1],w[2]+v[0]/c,w[3]+v[1]/c]);

	def pairingWith(self,circle2):
		w = self.w;
		w2 = circle2.w;
		Q = matrix(4,4,[0,-1,0,0, -1,0,0,0, 0,0,2,0, 0,0,0,2])/2;
		return w * Q * w2;

class CircleInSphere:
	def __init__(self,wPlus=None,m=None,alpha=None):
		'''
		Either spherical coordinates wPlus must be given, or
		a PointOnSphere m and an angle alpha.
		'''
		
		if m!=None and alpha!=None:
			self.m = m; #A PointOnSphere
			self.alpha = alpha;
			if wPlus==None:
				self.wPlus = vector([cot(alpha),m.z/sin(alpha),m.x/sin(alpha),m.y/sin(alpha)]);
			else:
				self.wPlus = wPlus;
		elif wPlus!=None:
			self.wPlus = wPlus;
			if alpha==None:
				self.alpha = arccot(wPlus[0]);
			else:
				self.alpha = alpha;
			if m==None:
				sinAlpha = sin(self.alpha);
				self.m = PointOnSphere(vector([sinAlpha*wPlus[2],sinAlpha*wPlus[3],sinAlpha*wPlus[1]]));
			else:
				self.m = m;
			
		else:
			print "Error: Not enough parameters to determine the given circle!";
			return;			

	def __str__(self):
		return "Circle: m = "+str(self.m)+", alpha = "+str(self.alpha)+".";

	def negative(self):
		wp = self.wPlus;
		neg_wp = vector([wp[0],-wp[1],-wp[2],-wp[3]]);
		#return CircleInSphere(wPlus=neg_wp,m=m.negative(),alpha=alpha);
		return CircleInSphere(wPlus=neg_wp);

	def toPlane(self):
		#m = self.m;
		#alpha = self.alpha;
		#X = m.x/(m.z+cos(alpha));
		#Y = m.y/(m.z+cos(alpha));
		#radius = sin(alpha)/(m.z+cos(alpha));
		#return CircleInPlane(PointOnPlane(vector([X,Y])),radius);

		wp = self.wPlus;
		w = vector([wp[0]-wp[1],wp[0]+wp[1],wp[2],wp[3]]);
		return CircleInPlane(w=w);

def numericMinimalPolynomial(a,maxDeg=16):
	for deg in range(1,maxDeg+1):
		minPoly = pslq([a^i for i in range(deg+1)]);
		if minPoly == None:
			continue;
		#Minimal polynomial found (at least approximately):
		R.<x> = QQ[];
		return sum([x^i*minPoly[i] for i in range(len(minPoly))]);
	#No minimal polynomial found up to the above degree:
	return "No minimal polynomial with degree <=",maxDeg,"found.";

########################################################################
### Tests: #############################################################
########################################################################

def test1():
	'''
	    2
	
	    3
	 0     1
	'''
	graphEmbedding = {0:[1,3,2],1:[2,3,0],2:[0,3,1],3:[0,1,2]};
	#radiiOfFixedNodes = {0:1, 1:1, 2:1};
	radiiOfFixedNodes = randomDict([0,1,2]);

	G = Graph(graphEmbedding);
	G.set_embedding(graphEmbedding);
	#G = Graph(4);
	#G.add_cycle([0,1,2]);
	#G.add_edge(0,3);
	#G.add_edge(1,3);
	#G.add_edge(2,3);
	#G.show();
	cp = CirclePacking(G,radiiOfFixedNodes);
	cp.show();
	cp.relationsAmongCurvatures();

def test2():
	'''
	      2
	 
	     4 3
	      5
	 0         1
	'''
	graphEmbedding = {0:[1,5,4,2],1:[2,3,5,0],2:[0,4,3,1],3:[2,4,5,1],4:[0,5,3,2],5:[1,3,4,0]};
	#radiiOfFixedNodes = {0:1, 1:2, 2:3};
	radiiOfFixedNodes = randomDict([0,1,2]);

	G = Graph(graphEmbedding);
	G.set_embedding(graphEmbedding);
	cp = CirclePacking(G,radiiOfFixedNodes);
	cp.show();
	cp.relationsAmongCurvatures();

def test3():
	'''
	Bipyramid over 5-gon: {2,6} * {0,1,5,4,3}
	           2
	
	           4
	      3         5
	           6
	 0                  1
	'''
	graphEmbedding = {0:[1,6,3,2],1:[5,6,0,2],2:[0,3,4,5,1],3:[0,6,4,2],4:[3,6,5,2],5:[4,6,1,2],6:[0,1,5,4,3]};
	#radiiOfFixedNodes = {0:1, 1:1, 2:1};
	radiiOfFixedNodes = randomDict([0,1,2]);

	G = Graph(graphEmbedding);
	#G.set_embedding(graphEmbedding);
	print "G is planar =",G.is_planar(set_embedding = True);
	cp = CirclePacking(G,radiiOfFixedNodes);
	cp.show();
	cp.relationsAmongCurvatures();

def test4():
	#G = graphs.CompleteGraph(4);
	#G = graphs.OctahedralGraph();
	G = graphs.IcosahedralGraph();
	#G = graphs.CycleGraph(5).join(Graph(2));
	#G.show();
	cp = CirclePacking(G);
	cp.show();
	cp.relationsAmongCurvatures();

def test5():
	triangulations = graphs.triangulations(4,minimum_connectivity=4);
	for G in triangulations:
		G.show();

def test6():
	#G = graphs.CompleteGraph(4);
	G = graphs.OctahedralGraph();
	#G = graphs.IcosahedralGraph();
	#G = graphs.CycleGraph(5).join(Graph(2));
	G.show();
	#findLinearRelationsAmongCurvatures(G);
	return findQuadraticRelationsAmongCurvatures(G);

#I = test6();

def test7():
	T = graphs.triangulations(8,minimum_connectivity=4);
	#T = [graphs.OctahedralGraph()];
	#T = [graphs.IcosahedralGraph()];
	for G in T:
		G.show();
		cp = CirclePacking(G);
		cp.show();
		findLinearRelationsAmongCurvatures(G);

#test7();
	
'''
OUTPUT:
Linear relations gibts bei (minimum_connectivity=4):
6_1 zwei Stueck (Octaeder)
8_2 4 
9_1 zwei Stueck (Zylinder ueber Dreieck, stacked an allen 3 Quadraten)
9_2 eine (
10_1 zwei
10_6 zwei
10_7 6
10_10 6
'''

def test8():
	#Variety given by circle packings:
	#G = graphs.TetrahedralGraph();
	G = graphs.OctahedralGraph();
	G.show();
	cpv = CirclePackingVariety(G);
	#cpv.checkLinearRelation(0,5,1,4);
	return cpv.A;

#A = test8();

def test09():
	G = graphs.OctahedralGraph();
	n = G.num_verts();
	varNames=None;
	#whichValues = lambda cp: [cp.curvatures[u] for u in cp.G.vertices()];
	#whichValues = lambda cp: [cp.R[u] for u in cp.G.vertices()];
	#whichValues = lambda cp: [cp.midpoints[u][0] for u in cp.G.vertices()] + [cp.midpoints[u][1] for u in cp.G.vertices()];

	whichValues = lambda cp: [cp.curvatures[u] for u in cp.G.vertices()] + [cp.midpoints[u][0] for u in cp.G.vertices()] + [cp.midpoints[u][1] for u in cp.G.vertices()];
	varNames = ['c'+str(i) for i in range(n)] + ['x'+str(i) for i in range(n)] + ['y'+str(i) for i in range(n)];

	#degrees = [1];
	degrees = [0,1,2];
	nSamples = 150;

	CirclePacking(G).show();
	return findRelationsInCirclePacking(G,degrees=degrees,whichValues=whichValues,nSamples=nSamples,varNames=varNames);

#I = test09();

def test09b():
	#G = graphs.OctahedralGraph();
	G = graphs.CycleGraph(6).join(Graph(2));
	print "G.vertices():",G.vertices();
	n = G.num_verts();
	varNames=None;
	whichValues = lambda cp: [cp.curvatures[u] for u in cp.G.vertices()];
	#whichValues = lambda cp: [cp.R[u] for u in cp.G.vertices()];
	#whichValues = lambda cp: [cp.midpoints[u][0] for u in cp.G.vertices()] + [cp.midpoints[u][1] for u in cp.G.vertices()];

	#whichValues = lambda cp: [cp.curvatures[u] for u in cp.G.vertices()] + [cp.midpoints[u][0] for u in cp.G.vertices()] + [cp.midpoints[u][1] for u in cp.G.vertices()];
	#arNames = ['c'+str(i) for i in range(n)] + ['x'+str(i) for i in range(n)] + ['y'+str(i) for i in range(n)];

	#degrees = [1];
	degrees = [0,1];
	nSamples = 15;

	G.show();
	CirclePacking(G).show();
	return findRelationsInCirclePacking(G,degrees=degrees,whichValues=whichValues,nSamples=nSamples,varNames=varNames);

#test09b();

def test10():
	R.<X,Y,w11,w12,w21,w22> = PolynomialRing(QQ,['X','Y','w11','w12','w21','w22']);
	print R;
	M = PointOnPlane(vector([X,Y]));
	print "M =",M;
	print "negM =",M.negative();
	print "invM =",M.inversion();
	
#test10();

def test11():
	def doTest(G):
		G.show()
		cp = CirclePacking(G);
		#cp.show()
		#print cp.gramMatrix;
		for x in cp.gramMatrix.column(0):
			#print x,":",pslq([1,x,x^2,x^3,x^4]);
			#print x,":",pslq([1,x,x^2]);
			print "Min poly for",x,":",numericMinimalPolynomial(x);

	#G = graphs.OctahedralGraph();
	#G = graphs.IcosahedralGraph();
	#doTest(G);
	#return;

	#T = [graphs.OctahedralGraph()];
	#T = [graphs.IcosahedralGraph()];

	T = graphs.triangulations(8,minimum_connectivity=4);
	for G in T:
		doTest(G);
	
#test11();
	

