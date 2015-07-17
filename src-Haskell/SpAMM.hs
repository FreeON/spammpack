module SpAMM
( MatrixTree
, treeAdd
, treeMult
, treeMultTol
, treeTranspose
) where

-- matrix algebra on MatrixTrees, including recursive multiplication that
-- returns Zeros when products of norms fall below tolerance (SpAMM)

import MatrixTree (addSubtreeNorms, ifZeroReplace, MatrixTree, MTree(..),
                   nextPowOf2, norm, setNorm, size, Size, valueNorm)

treeTranspose :: MatrixTree -> MatrixTree
treeTranspose (h, w, mTree) = (w, h, mTTrans mTree)

mTTrans :: MTree -> MTree
mTTrans (Square s x tl tr bl br) = Square s x (mTTrans tl) (mTTrans bl)
                                              (mTTrans tr) (mTTrans br)
mTTrans tree = tree

treeAdd :: MatrixTree -> MatrixTree -> MatrixTree
treeAdd (h1, w1, mTree1) (h2, w2, mTree2) =
        if h1 == h2 && w1 == w2 then (h1, w1, mTAdd mTree1 mTree2)
        else error "matrices don't match for addition"

mTAdd :: MTree -> MTree -> MTree
mTAdd (Zero _) mTree = mTree
mTAdd mTree (Zero _) = mTree
mTAdd (Leaf _ x) (Leaf _ y) =
      setNorm . ifZeroReplace $ Leaf 0 (x + y)
mTAdd (Square s _ tl1 tr1 bl1 br1) (Square _ _ tl2 tr2 bl2 br2) =
      ifZeroReplace $ Square s n tl tr bl br
      where tl = tl1 `mTAdd` tl2 ; tr = tr1 `mTAdd` tr2
            bl = bl1 `mTAdd` bl2 ; br = br1 `mTAdd` br2
            n = addSubtreeNorms [tl, tr, bl, br]
mTAdd _ _ = error "matrices don't match for addition"

treeMult :: MatrixTree -> MatrixTree -> MatrixTree
treeMult = treeMultTol 0

treeMultTol :: Double -> MatrixTree -> MatrixTree -> MatrixTree
treeMultTol tol (h1, w1, mTree1) (h2, w2, mTree2) =
            if w1 == h2 then cutToSize (h1, w2, mTMult tol expMTree1 expMTree2)
            else error "matrices don't match for multiplication"
            where expMTree1 = expandMTree mTree1 (size mTree2)
                  expMTree2 = expandMTree mTree2 (size mTree1)

mTMult :: Double -> MTree -> MTree -> MTree
mTMult _ zero@(Zero _) _ = zero
mTMult _ _ zero@(Zero _) = zero
mTMult tol (Leaf m x) (Leaf n y) = if m * n <= tol then Zero 1
                                   else setNorm $ Leaf 0 (x * y)
mTMult tol (Square s m tl1 tr1 bl1 br1) (Square _ n tl2 tr2 bl2 br2) =
         if m * n <= tol then Zero s
         else ifZeroReplace $ Square s p tl tr bl br
         where tl = (tl1 `mTTimes` tl2) `mTAdd` (tr1 `mTTimes` bl2)
               tr = (tl1 `mTTimes` tr2) `mTAdd` (tr1 `mTTimes` br2)
               bl = (bl1 `mTTimes` tl2) `mTAdd` (br1 `mTTimes` bl2)
               br = (bl1 `mTTimes` tr2) `mTAdd` (br1 `mTTimes` br2)
               mTTimes = mTMult tol
               p = addSubtreeNorms [tl, tr, bl, br]
mTMult _ _ _ = error "matrices don't match for multiplication"

expandMTree :: MTree -> Size -> MTree
expandMTree mTree n = if size mTree >= n then mTree
                      else expandMTree (Square dbl m mTree zro zro zro) n
                      where dbl = 2 * size mTree
                            m = norm mTree
                            zro = Zero (size mTree)

cutToSize :: MatrixTree -> MatrixTree
cutToSize tree@(h, w, mTree) = if size mTree <= nextPowOf2 (max h w) then tree
                               else cutToSize (h, w, cut mTree)
                               where cut (Zero s) = Zero (s `div` 2)
                                     cut (Square _ _ tl _ _ _) = tl