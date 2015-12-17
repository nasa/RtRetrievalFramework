# Routine for find the Longest Common Substring of two strings
import sys
import random

# kludge: infinity is a very large number
inf = sys.maxint

# Define a class for a node in the suffix tree
class SuffixNode(dict):
    def __init__(self):
        self.suffixLink = None # Suffix link as defined by Ukkonen

    def Print(self,str,ws=""):
        for t in self:
            k,p,s = self[t]
            if p == inf:
                print "%s%s" % (ws, str[k:])
            else:
                print "%s%s" % (ws, str[k:p+1])
                s.Print(str,ws+"|"*(p-k+1))
        

class LCS:
    def __init__(self,str1,str2):

        check_str = str1 + str2
        inf = len(check_str)
        self.check_str = check_str   #Keep a reference to check_str to ensure the string is not garbage collected
        self.seed = SuffixNode() #Seed is a dummy node. Suffix link of root points to seed. For any char,there is a link from seed to root
        self.root = SuffixNode() # Root of the suffix tree
        self.root.suffixLink = self.seed
        self.root.depth = 0
        self.deepest = 0,0

        # For each character of check_str[i], create suffixtree for check_str[0:i]
        s = self.root; k=0
        for i in range(len(check_str)):
            self.seed[check_str[i]] = -2,-2,self.root
            oldr = self.seed
            t = check_str[i]
            #Traverse the boundary path of the suffix tree for check_str[0:i-1]
            while True:
               # Decend the suffixtree until state s has a transition for the stringstr[k:i-1]
                while i>k:
                   kk,pp,ss = s[check_str[k]]
                   if pp-kk < i-k:
                       k = k + pp-kk+1
                       s = ss
                   else:
                       break
               # Exit this loop if s has a transition for the string check_str[k:i] (itmeans check_str[k:i] is repeated);
               # Otherwise, split the state if necessary
                if i>k:
                    tk = check_str[k]
                    kp,pp,sp = s[tk]
                    if t == check_str[kp+i-k]:
                        break
                    else: # Split the node
                        r = SuffixNode()
                        j = kp+i-k
                        tj = check_str[j]
                        r[tj] = j, pp, sp
                        s[check_str[kp]] = kp,j-1, r
                        r.depth = s.depth + (i-k)
                        sp.depth = r.depth + pp - j + 1
                        if j<len(str1)<i and r.depth>self.deepest[0]:
                            self.deepest = r.depth,j-1
                elif s.has_key(t):
                    break
                else:
                    r = s
               # Add a transition from r that starts with the letter check_str[i]
                tmp = SuffixNode()
                r[t] = i,inf,tmp
                # Prepare for next iteration
                oldr.suffixLink = r
                oldr = r
                s = s.suffixLink
            # Last remaining endcase
            oldr.suffixLink = s

    def LongestCommonSubstring(self):
        return self.check_str[self.deepest[1]-self.deepest[0]+1:self.deepest[1]+1]

    def Print(self):
        self.root.Print(self.check_str)
