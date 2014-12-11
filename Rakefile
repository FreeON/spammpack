require 'rubygems'
require 'rake'
require 'rdoc'
require 'date'
require 'yaml'
require 'tmpdir'
require 'jekyll'

desc "Generate blog files"
task :generate do
  Jekyll::Site.new(Jekyll.configuration({
    "source"      => ".",
    "destination" => "_site"
  })).process
end


desc "Generate and publish blog to gh-pages"
task :publish => [:generate] do
  Dir.mktmpdir do |tmp|
    system "mv _site/* #{tmp}"
    system "git remote update"
    system "git stash"
    system "git checkout gh-pages"
    system "git reset --hard origin/gh-pages"
    system "git rm -r ."
    system "rm -rf *"
    system "mv #{tmp}/* ."
    system "touch .nojekyll"
    system "git add .nojekyll"
    system "git add ."
    message = "Site updated at #{Time.now.utc}"
    system "git commit -am #{message.shellescape}"
    system "git push"
    system "git checkout master"
    system "git stash pop"
  end
end

task :default => :generate
