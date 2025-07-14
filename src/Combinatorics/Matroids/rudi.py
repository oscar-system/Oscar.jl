def U24minors(M):
        if M.full_rank()<2 or M.full_corank()<2:
            return
        for F in M.flats(M.full_rank()-2):
            H = []
            G = set(M.groundset()-F)
            while G:
                e = G.pop()
                H.append(e)
                G = G-M.closure(F|set([e]))
            
            if len(H)<4:
                continue
            
            I = M.max_independent(F)
            for i in range(len(H)):
                a=H[i]
                for j in range(i):
                    b=H[j]
                    for k in range(j):
                        c=H[k]
                        for l in range(k):
                            d=H[l]
                            yield I,a,b,c,d

class TutteGroup:
    def __init__(self, M, char2 = False):
        
        B = list(M.bases())
        idx = {B[i] : i for i in range(len(B))}
        idx[frozenset()] = len(B)
        self._M = M
        self._bases = B
        self._idx = idx
        
        if char2:
            R = [{frozenset():1}]
        else:
            R = [{frozenset():2}]
        
        
        for X in M.nonbases():
            if M.rank(X)==M.full_rank()-1: 
                C=set(M.circuit(X))
                D=set(M.cocircuit(M.groundset()-X))
                e = C.pop()
                f = D.pop()
                for g in C:
                    for h in D:
                        R.append(self.cross_ratio(X-set([e,g]),e,g,f,h))
        self._TT = matrix(ZZ, 1,len(B)+1, sparse = True)
        
        self.add_relations(R)
        
    def add_relations(self, R):     
        T  = matrix(ZZ, len(R)+1,len(self._bases)+1, sparse = True)
        for i in range(len(R)):
            for X, e in R[i].items():
                T[i, self._idx[X]]=e  
        T = self._TT.stack(T).sparse_matrix()
        TT = T.hermite_form()
        self._piv = TT.pivots()
        self._tr = len(self._piv)
        self._TT = TT.matrix_from_rows_and_columns(range(self._tr), range(TT.ncols())).stack(vector(ZZ,TT.ncols())).sparse_matrix()
        self._free = frozenset([self._bases[i] for i in range(self._TT.ncols()-1) if self._TT.column(i).is_zero()])
    
    def __repr__(self):
        return 'Tutte group of ' + repr(self._M)
      
    def cross_ratio(self,I,a,b,c,d):
        return {I|set([a,c]):1,
                I|set([b,d]):1,
                I|set([a,d]):-1,
                I|set([b,c]):-1,
                frozenset():sum([a<c, c<b, b<d, d<a])}
    
    def is_zero(self, t):
        for B in t.keys():
            if B in self._free:
                return False
        r = self._tr
        for X, e in t.items():
            self._TT[r, self._idx[X]]=e
        i = 0
        jj = 0
        for j in self._piv:
            if self._TT[r, j]:
                if self._TT[r,j]%self._TT[i,j]:
                    self._TT.set_row_to_multiple_of_row(r,r,0)
                    return False
                self._TT.add_multiple_of_row(r, i, -self._TT[r,j]//self._TT[i,j])   
            i += 1
        if not self._TT.row(r).is_zero():
            self._TT.set_row_to_multiple_of_row(r,r,0)
            return False
        return True 
    
    def cross_ratios(self):
        for I,a,b,c,d in U24minors(self._M):
            yield self.cross_ratio(I,a,b,c,d)
            yield self.cross_ratio(I,a,c,d,b)
            yield self.cross_ratio(I,a,d,b,c)
            
    
    def DW_condition(self):
        for cr in self.cross_ratios():
            if self.is_zero(cr):
                return False
        return True
    
    def is_char2(self):
        return self.is_zero({frozenset():1})
    
    def ellipse(self,I,a,b,c,d,e,f):
        R = self.cross_ratio(I|set([e]),a,b,c,d) 
        R2 = self.cross_ratio(I|set([f]),a,b,c,d)
        for X, e in R2.items():
            if X not in R:
                R[X]=-e
            else:
                R[X]-=e
        return R
            
    def contract(self, e):
        B = [b for b in self._bases if e in b]
        B2 = [b for b in self._bases if e not in b]
        T = self._TT.matrix_from_rows_and_columns(range(self._TT.nrows()), [self._idx[b] for b in B2+B+[frozenset()]]).hermite_form()
        T = T.matrix_from_rows_and_columns([r for r in range(len(T.pivots())) if T.pivots[r]>=len(B2)], range(len(B2), len(B2)+len(B)+1))                                  
        idx = { B[i]:i for i in range(len(B))}
        idx[frozenset()]=len(B)
        T = T.hermite_form()
        self._M = self._M/e
        self._bases = B
        self._idx = idx
        self._TT = T.stack(matrix(ZZ, 1,len(B)+1, sparse = True))
        self._piv = T.pivots()
        self._tr = len(self._piv)
        self._free = frozenset([self._bases[i] for i in range(T.ncols()-1) if T.column(i).is_zero()])
    
    def delete(self, e):
        B = [b for b in self._bases if e not in b]
        B2 = [b for b in self._bases if e in b]
        T = self._TT.matrix_from_rows_and_columns(range(self._TT.nrows()), [self._idx[b] for b in B2+B+[frozenset()]]).hermite_form()
        T = T.matrix_from_rows_and_columns([r for r in range(len(T.pivots())) if T.pivots()[r]>=len(B2)], range(len(B2), len(B2)+len(B)+1))                                  
        idx = { B[i]:i for i in range(len(B))}
        idx[frozenset()]=len(B)
        T = T.hermite_form()
        self._M = self._M\e
        self._bases = B
        self._idx = idx
        self._TT = T.stack(matrix(ZZ, 1,len(B)+1, sparse = True))
        self._piv = T.pivots()
        self._tr = len(self._piv)
        self._free = frozenset([self._bases[i] for i in range(T.ncols()-1) if T.column(i).is_zero()])
