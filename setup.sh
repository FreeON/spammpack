#!/bin/bash

echo "Installing bundler and all necessary gems..."
gem install bundler
bundle install --binstubs

echo "The gems can be updated by running"
echo "bundle update"
