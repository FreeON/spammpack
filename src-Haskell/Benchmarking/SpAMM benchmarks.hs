import Control.DeepSeq (NFData(..))
import Criterion.Main
import MatrixTree
import SpAMM
import System.IO (hClose, hPutStrLn, openFile, IOMode(WriteMode))
import System.Process

instance NFData MatrixTree

main = do let sizes = fmap (2^) [11..20]
          let fstFilePaths = fmap (\n -> "fst" ++ show n ++ ".txt") sizes
          let sndFilePaths = fmap (\n -> "snd" ++ show n ++ ".txt") sizes
          mapM_ (uncurry createMatrix) (zip fstFilePaths sizes)
          mapM_ (uncurry createMatrix) (zip sndFilePaths sizes)
          fstTrees <- mapM readTreeFromMatrixMarket fstFilePaths
          sndTrees <- mapM readTreeFromMatrixMarket sndFilePaths
          callProcess "rm" $ fstFilePaths ++ sndFilePaths
          let sizeTreeGroups = zip sizes $ zip fstTrees sndTrees
          defaultMain [bgroup "add"  $ fmap (makeBench treeAdd)  sizeTreeGroups,
                       bgroup "mult" $ fmap (makeBench treeMult) sizeTreeGroups]

createMatrix :: FilePath -> Int -> IO ()
createMatrix filePath size = do
             handle <- openFile filePath WriteMode
             matrix <- readProcess "python" ["../generate-matrix.py", "-N " ++ show size] []
             hPutStrLn handle matrix
             hClose handle

makeBench :: (MatrixTree -> MatrixTree -> MatrixTree) -> (Int, (MatrixTree, MatrixTree)) -> Benchmark
makeBench op (size, (tree1, tree2)) = bench (show size) $ nf (op tree1) tree2