module SpAMM
( MatrixPair
, pairAdd
, pairMult
, pairMultTol
, pairTranspose
, treeAdd
, treeMult
, treeTranspose
) where

-- matrix algebra on MatrixPairs, including matrix multiplication that recurs on subtrees
-- and returns Zeros when products of norms fall below tolerance (SpAMM)

import MatrixPair (addSubtreeNorms, height, ifZeroReplace, MatrixPair, MatrixTree(..),
                   norm, valueNorm, width)

pairTranspose :: MatrixPair -> MatrixPair
pairTranspose (i, j, tree) = (j, i, treeTranspose tree)

treeTranspose :: MatrixTree -> MatrixTree
treeTranspose (Zero h w)               = Zero w h
treeTranspose tree@(Value _ _)         = tree
treeTranspose (Rect h w x tl tr bl br) = Rect w h x
                                              (treeTranspose tl) (treeTranspose bl)
                                              (treeTranspose tr) (treeTranspose br)

pairAdd :: MatrixPair -> MatrixPair -> MatrixPair
pairAdd (h1, w1, tree1) (h2, w2, tree2) = if h1 == h2 && w1 == w2
                                          then (h1, w1, treeAdd tree1 tree2)
                                          else error "matrices don't match for addition"

treeAdd :: MatrixTree -> MatrixTree -> MatrixTree

treeAdd (Zero _ _) tree = tree

treeAdd tree (Zero _ _) = tree

treeAdd (Value _ x) (Value _ y) = if x + y == 0 then Zero 1 1
                                  else Value (valueNorm $ x + y) (x + y)

treeAdd (Rect h w _ tl1 tr1 bl1 br1) (Rect _ _ _ tl2 tr2 bl2 br2) =
         ifZeroReplace $ Rect h w x tlsum trsum blsum brsum
         where tlsum = tl1 `treeAdd` tl2 ; trsum = tr1 `treeAdd` tr2
               blsum = bl1 `treeAdd` bl2 ; brsum = br1 `treeAdd` br2
               x = addSubtreeNorms . fmap norm $ [tlsum, trsum, blsum, brsum]

treeAdd _ _ = error "matrices don't match for addition"

pairMult :: MatrixPair -> MatrixPair -> MatrixPair
pair1 `pairMult` pair2 = pairMultTol pair1 pair2 0

pairMultTol :: MatrixPair -> MatrixPair -> Double -> MatrixPair
pairMultTol (h1, w1, tree1) (h2, w2, tree2) tol =
            if w1 == h2 then (h1, w2, treeMult tol tree1 tree2)
            else error "matrices don't match for multiplication"

treeMult :: Double -> MatrixTree -> MatrixTree -> MatrixTree

treeMult _ (Zero h _) tree = Zero h (width tree)

treeMult _ tree (Zero _ w) = Zero (height tree) w

treeMult tol tree1@(Rect _ _ _ _ _ _ _) tree2@(Value _ _) =
             treeTranspose $ treeMult tol (treeTranspose tree2) (treeTranspose tree1)

treeMult tol (Value m x) (Value n y) = if m * n <= tol then Zero 1 1
                                       else Value (valueNorm $ x * y) (x * y)

treeMult tol val@(Value m _) (Rect _ w n tl tr bl br) =
         if m * n <= tol then Zero 1 w
         else ifZeroReplace $ Rect 1 w x tlmult trmult blmult brmult
         where tlmult = val `treeTimes` tl ; trmult = val `treeTimes` tr
               blmult = val `treeTimes` bl ; brmult = val `treeTimes` br
               treeTimes = treeMult tol
               x = addSubtreeNorms . fmap norm $ [tlmult, trmult, blmult, brmult]

treeMult tol (Rect h1 _ m tl1 tr1 bl1 br1) (Rect _ w2 n tl2 tr2 bl2 br2) =
         if m * n <= tol then Zero h1 w2
         else ifZeroReplace $ Rect h1 w2 x tlmult trmult blmult brmult
         where tlmult = (tl1 `treeTimes` tl2) `treeAdd` (tr1 `treeTimes` bl2)
               trmult = (tl1 `treeTimes` tr2) `treeAdd` (tr1 `treeTimes` br2)
               blmult = (bl1 `treeTimes` tl2) `treeAdd` (br1 `treeTimes` bl2)
               brmult = (bl1 `treeTimes` tr2) `treeAdd` (br1 `treeTimes` br2)
               treeTimes = treeMult tol
               x = addSubtreeNorms . fmap norm $ [tlmult, trmult, blmult, brmult]
