# Installation for Ubuntu

## Install jekyll
sudo apt-get install ruby ruby-dev build-essential
## Update paths
echo '# Install Ruby Gems to ~/gems' >> ~/.bashrc
echo 'export GEM_HOME=$HOME/gems' >> ~/.bashrc
echo 'export PATH=$HOME/gems/bin:$PATH' >> ~/.bashrc
source ~/.bashrc
## Install Jekyll and bundler
gem install jekyll bundler

# Getting Updates 
bundle init
bundle install
bundle update

