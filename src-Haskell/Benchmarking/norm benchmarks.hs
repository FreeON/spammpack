import Control.DeepSeq (NFData(..))
import Criterion.Main
import MatrixTree
import SpAMM
import System.Random (getStdGen, randomRs, StdGen)

instance NFData MTree

main = do let sizes = fmap (2^) [5..11]
          gen <- getStdGen
          let trees = fmap (makeRandomMTree gen) sizes
          let sizeTreePairs = zip sizes trees
          defaultMain [bgroup "norm" $ fmap makeBench sizeTreePairs]

makeRandomMTree :: StdGen -> Int -> MTree
makeRandomMTree gen size = mTree
                where (_, _, mTree) = matrixListToTree (size, size, ijxs)
                      ijxs = zipWith (\(i, j) x -> (i, j, x)) indices randomNums
                      indices = [(i, j) | j <- [1..size], i <- [1..size]]
                      randomNums = take (size^2) $ randomRs (0.0, 1.0) gen

makeBench :: (Int, MTree) -> Benchmark
makeBench (size, mTree) = bench (show size) $ nf setNorm mTree