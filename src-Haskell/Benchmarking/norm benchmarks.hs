import Control.DeepSeq (NFData(..))
import Criterion.Main
import MatrixTree
import SpAMM
import System.IO (hClose, hPutStrLn, openFile, IOMode(WriteMode))
import System.Process

instance NFData MatrixTree

main = do let sizes = fmap (2^) [1..10]
          let filePaths = fmap (\n -> show n ++ ".txt") sizes
          mapM_ (uncurry createMatrix) (zip filePaths sizes)
          trees <- mapM readTreeFromMatrixMarket filePaths
          callProcess "rm" filePaths
          let sizeTreePairs = zip sizes trees
          defaultMain [bgroup "norm" $ fmap makeBench sizeTreePairs]

createMatrix :: FilePath -> Int -> IO ()
createMatrix filePath size = do
             handle <- openFile filePath WriteMode
             matrix <- readProcess "python" ["../generate-matrix.py", "-N " ++ show size] []
             hPutStrLn handle matrix
             hClose handle

makeBench :: (Int, MatrixTree) -> Benchmark
makeBench (size, tree) = bench (show size) $ nf setNorm tree