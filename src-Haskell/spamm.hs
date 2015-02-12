module SpAMM
( MatrixTree
, treeAdd
, treeMult
, treeMultTol
, treeTranspose
) where

-- matrix algebra on MatrixTrees, including matrix multiplication that recurs on subtrees
-- and returns Zeros when products of norms fall below tolerance (SpAMM)

import MatrixTree (addSubtreeNorms, getNorm, ifZeroReplace, MatrixTree(..), rectOrder, valueNorm)

treeTranspose :: MatrixTree -> MatrixTree
treeTranspose (Zero t l h w)               = Zero l t w h
treeTranspose (Value i j x y)              = Value j i x y
treeTranspose (Rect t l h w x tl tr bl br) = Rect l t w h x
                                                  (treeTranspose tl) (treeTranspose bl)
                                                  (treeTranspose tr) (treeTranspose br)

treeAdd :: MatrixTree -> MatrixTree -> MatrixTree

treeAdd (Zero i j k l) tree@(Zero m n p q)
        | [i,j,k,l] == [m,n,p,q] = tree
        | otherwise              = error "matrices don't match for addition"

treeAdd tree zeroTree@(Zero _ _ _ _) = treeAdd zeroTree tree

treeAdd (Zero t l h w) tree@(Value i j _ _)
        | [t,l,h,w] == [i,j,1,1] = tree
        | otherwise              = error "matrices don't match for addition"

treeAdd (Zero t l h w) tree@(Rect u m j v _ _ _ _ _)
        | [t,l,h,w] == [u,m,j,v] = tree
        | otherwise              = error "matrices don't match for addition"

treeAdd (Value i j _ x) (Value k l _ y)
        | [i,j] == [k,l] = if x + y == 0 then Zero i j 1 1
                           else Value i j (valueNorm $ x + y) (x + y)
        | otherwise      = error "matrices don't match for addition"

treeAdd rect1@(Rect t l h w _ _ _ _ _) rect2@(Rect u m j v _ _ _ _ _)
        | [t,l,h,w] == [u,m,j,v] = ifZeroReplace (Rect t l h w x tlsum trsum blsum brsum)
        | otherwise              = error "matrices don't match for addition"
        where Rect _ _ _ _ _ tl1 tr1 bl1 br1 = rectOrder rect1
              Rect _ _ _ _ _ tl2 tr2 bl2 br2 = rectOrder rect2
              tlsum = tl1 `treeAdd` tl2 ; trsum = tr1 `treeAdd` tr2
              blsum = bl1 `treeAdd` bl2 ; brsum = br1 `treeAdd` br2
              x = addSubtreeNorms . fmap getNorm $ [tlsum, trsum, blsum, brsum]

treeAdd _ _ = error "matrices don't match for addition"

treeMult :: MatrixTree -> MatrixTree -> MatrixTree
tree1 `treeMult` tree2 = treeMultTol tree1 tree2 0

treeMultTol :: MatrixTree -> MatrixTree -> Double -> MatrixTree

treeMultTol (Zero t1 l1 h1 w1) (Zero t2 l2 h2 w2) _
            | [l1,w1] == [t2,h2] = Zero t1 l2 h1 w2
            | otherwise          = error "matrices don't match for multiplication"

treeMultTol tree zeroTree@(Zero _ _ _ _) _ =
            treeTranspose $ treeMult (treeTranspose zeroTree) (treeTranspose tree)

treeMultTol (Zero t l h w) (Value i j _ _) _
            | l == i && w `elem` [0,1] = Zero t j h w
            | otherwise                = error "matrices don't match for multiplication"

treeMultTol (Zero t1 l1 h1 w1) (Rect t2 l2 h2 w2 _ _ _ _ _) _
            | [l1,w1] == [t2,h2] = Zero t1 l2 h1 w2
            | otherwise          = error "matrices don't match for multiplication"

treeMultTol tree1@(Rect _ _ _ _ _ _ _ _ _) tree2@(Value _ _ _ _) tol =
            treeTranspose $ treeMultTol (treeTranspose tree2) (treeTranspose tree1) tol

treeMultTol (Value i j m x) (Value k l n y) tol
            | j == k    = if m * n <= tol then Zero i l 1 1
                          else Value i l (valueNorm $ x * y) (x * y)
            | otherwise = error "matrices don't match for multiplication"

treeMultTol val@(Value i j m _) (Rect t l h w n tl tr bl br) tol
            | [j,1] == [t,h] = if m * n <= tol then Zero i l 1 w
                               else ifZeroReplace (Rect i l 1 w x
                                                        tlmult trmult blmult brmult)
            | otherwise      = error "matrices don't match for multiplication"
            where [tlmult, trmult, blmult, brmult] = fmap multiply [tl, tr, bl, br]
                  multiply tree = treeMultTol val tree tol
                  x = addSubtreeNorms . fmap getNorm $ [tlmult, trmult, blmult, brmult]

treeMultTol (Rect t1 l1 h1 w1 m tl1 tr1 bl1 br1) (Rect t2 l2 h2 w2 n tl2 tr2 bl2 br2) tol
            | [l1,w1] == [t2,h2] = if m * n <= tol then Zero t1 l2 h1 w2
                                   else ifZeroReplace (Rect t1 l2 h1 w2 x
                                                            tlmult trmult blmult brmult)
            | otherwise          = error "matrices don't match for multiplication"
            where tlmult = treeAdd (treeMultTol tl1 tl2 tol) (treeMultTol tr1 bl2 tol)
                  trmult = treeAdd (treeMultTol tl1 tr2 tol) (treeMultTol tr1 br2 tol)
                  blmult = treeAdd (treeMultTol bl1 tl2 tol) (treeMultTol br1 bl2 tol)
                  brmult = treeAdd (treeMultTol bl1 tr2 tol) (treeMultTol br1 br2 tol)
                  x = addSubtreeNorms . fmap getNorm $ [tlmult, trmult, blmult, brmult]