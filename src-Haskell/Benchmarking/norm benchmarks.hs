import Control.DeepSeq (NFData(..))
import Criterion.Main
import MatrixTree
import SpAMM
import System.Random (getStdGen, randomRs, StdGen)

instance NFData MatrixTree

main = do let sizes = fmap (2^) [5..11]
          gen <- getStdGen
          let trees = fmap (flip makeRandomTree gen) sizes
          let sizeTreePairs = zip sizes trees
          defaultMain [bgroup "norm" $ fmap makeBench sizeTreePairs]

makeRandomTree :: Int -> StdGen -> MatrixTree
makeRandomTree size gen = indexedListToTree (size, size, ijxs)
               where ijxs = zipWith (\(i, j) x -> (i, j, x)) indices randomNums
                     indices = [(i, j) | j <- [1..size], i <- [1..size]]
                     randomNums = take (size^2) $ randomRs (0.0, 1.0) gen

makeBench :: (Int, MatrixTree) -> Benchmark
makeBench (size, tree) = bench (show size) $ nf setNorm tree