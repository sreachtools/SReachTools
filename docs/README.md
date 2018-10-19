# Docs for SReachTools

We use Github pages and Jekyll to power the documentation website for [SReachTools](https://unm-hscl.github.io/SReachTools/).

## Installation for Ubuntu

1. Install jekyll 
    ```
    sudo apt-get install ruby ruby-dev build-essential
    ```
1. Update paths 
    ```
    echo '# Install Ruby Gems to ~/gems' >> ~/.bashrc
    echo 'export GEM_HOME=$HOME/gems' >> ~/.bashrc
    echo 'export PATH=$HOME/gems/bin:$PATH' >> ~/.bashrc
    source ~/.bashrc
    ```
1. Install Jekyll and bundler
    ```
    gem install jekyll bundler
    ```
1. Getting Updates 
    ```
    bundle init
    bundle install
    bundle update
    ```
