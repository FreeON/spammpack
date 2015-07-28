import Control.DeepSeq (NFData(..))
import Criterion.Main
import MatrixTree
import SpAMM
import System.Random (getStdGen, newStdGen, randomRs, StdGen)

instance NFData MTree

main = do let sizes = fmap (2^) [5..10]
          gen <- getStdGen
          let fstTrees = fmap (flip makeRandomTree gen) sizes
          gen <- newStdGen
          let sndTrees = fmap (flip makeRandomTree gen) sizes
          let sizeTreeGroups = zip sizes $ zip fstTrees sndTrees
          defaultMain [bgroup "add"  $ fmap (makeBench treeAdd)  sizeTreeGroups,
                       bgroup "mult" $ fmap (makeBench treeMult) sizeTreeGroups]

makeRandomTree :: Int -> StdGen -> MatrixTree
makeRandomTree size gen = matrixListToTree (size, size, ijxs)
               where ijxs = zipWith (\(i, j) x -> (i, j, x)) indices randomNums
                     indices = [(i, j) | j <- [1..size], i <- [1..size]]
                     randomNums = take (size^2) $ randomRs (0.0, 1.0) gen

makeBench :: (MatrixTree -> MatrixTree -> MatrixTree) -> (Int, (MatrixTree, MatrixTree)) -> Benchmark
makeBench op (size, (tree1, tree2)) = bench (show size) $ nf (op tree1) tree2