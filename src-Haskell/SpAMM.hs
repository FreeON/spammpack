module SpAMM
( MatrixTree
, treeAdd
, treeMult
, treeMultTol
, treeTranspose
) where

-- matrix algebra on MatrixTrees, including matrix multiplication that recurs on subtrees
-- and returns Zeros when products of norms fall below tolerance (SpAMM)

import MatrixTree (addSubtreeNorms, height, ifZeroReplace, MatrixTree(..),
                   norm, rectOrder, valueNorm, width)

treeTranspose :: MatrixTree -> MatrixTree
treeTranspose (Zero t l h w)               = Zero l t w h
treeTranspose (Value i j x y)              = Value j i x y
treeTranspose (Rect t l h w x tl tr bl br) = Rect l t w h x
                                                  (treeTranspose tl) (treeTranspose bl)
                                                  (treeTranspose tr) (treeTranspose br)

treeAdd :: MatrixTree -> MatrixTree -> MatrixTree
treeAdd tree1 tree2 = if tree1 `matchForAdd` tree2 then treePlus tree1 tree2
                      else error "matrices don't match for addition"

matchForAdd :: MatrixTree -> MatrixTree -> Bool
matchForAdd tree1 tree2 = top tree1    == top tree2    && left tree1  == left tree2  &&
                          height tree1 == height tree2 && width tree1 == width tree2

treePlus :: MatrixTree -> MatrixTree -> MatrixTree

treePlus (Zero _ _ _ _) tree = tree

treePlus tree (Zero _ _ _ _) = tree

treePlus (Value i j _ x) (Value _ _ _ y) = if x + y == 0 then Zero i j 1 1
                                           else Value i j (valueNorm $ x + y) (x + y)

treePlus rect1@(Rect t l h w _ _ _ _ _) rect2 =
         ifZeroReplace (Rect t l h w x tlsum trsum blsum brsum)
         where Rect _ _ _ _ _ tl1 tr1 bl1 br1 = rectOrder rect1
               Rect _ _ _ _ _ tl2 tr2 bl2 br2 = rectOrder rect2
               tlsum = tl1 `treePlus` tl2 ; trsum = tr1 `treePlus` tr2
               blsum = bl1 `treePlus` bl2 ; brsum = br1 `treePlus` br2
               x = addSubtreeNorms . fmap norm $ [tlsum, trsum, blsum, brsum]

treePlus _ _ = error "matrices don't match for addition"

treeMult :: MatrixTree -> MatrixTree -> MatrixTree
tree1 `treeMult` tree2 = treeMultTol tree1 tree2 0

treeMultTol :: MatrixTree -> MatrixTree -> Double -> MatrixTree
treeMultTol tree1 tree2 tol = if tree1 `matchForMult` tree2 then treeTimesTol tree1 tree2 tol
                              else error "matrices don't match for multiplication"

matchForMult :: MatrixTree -> MatrixTree -> Bool
matchForMult tree1 tree2 = left tree1 == top tree2 && width tree1 == height tree2

treeTimesTol :: MatrixTree -> MatrixTree -> Double -> MatrixTree

treeTimesTol (Zero t _ h w) tree _ = Zero t (left tree) h (min w $ width tree)

treeTimesTol tree (Zero _ l h w) _ = Zero (top tree) l (min h $ height tree) w

treeTimesTol tree1@(Rect _ _ _ _ _ _ _ _ _) tree2@(Value _ _ _ _) tol =
             treeTranspose $ treeTimesTol (treeTranspose tree2) (treeTranspose tree1) tol

treeTimesTol (Value i _ m x) (Value _ l n y) tol = if m * n <= tol then Zero i l 1 1
                                                   else Value i l (valueNorm $ x * y) (x * y)

treeTimesTol val@(Value i _ m _) (Rect _ l _ w n tl tr bl br) tol =
             if m * n <= tol then Zero i l 1 w
             else ifZeroReplace (Rect i l 1 w x tlmult trmult blmult brmult)
             where tlmult = treeTimesTol val tl tol ; trmult = treeTimesTol val tr tol
                   blmult = treeTimesTol val bl tol ; brmult = treeTimesTol val br tol
                   x = addSubtreeNorms . fmap norm $ [tlmult, trmult, blmult, brmult]

treeTimesTol (Rect t1 _ h1 _ m tl1 tr1 bl1 br1) (Rect _ l2 _ w2 n tl2 tr2 bl2 br2) tol =
             if m * n <= tol then Zero t1 l2 h1 w2
             else ifZeroReplace (Rect t1 l2 h1 w2 x tlmult trmult blmult brmult)
             where tlmult = treePlus (treeTimesTol tl1 tl2 tol) (treeTimesTol tr1 bl2 tol)
                   trmult = treePlus (treeTimesTol tl1 tr2 tol) (treeTimesTol tr1 br2 tol)
                   blmult = treePlus (treeTimesTol bl1 tl2 tol) (treeTimesTol br1 bl2 tol)
                   brmult = treePlus (treeTimesTol bl1 tr2 tol) (treeTimesTol br1 br2 tol)
                   x = addSubtreeNorms . fmap norm $ [tlmult, trmult, blmult, brmult]
