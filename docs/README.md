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
    bundle update
    bundle init
    bundle install
    ```
    In case, you receive a cryptic error of the form:
```
Traceback (most recent call last):
	2: from /home/abyvinod/gems/bin/bundler:23:in `<main>'
	1: from /usr/lib/ruby/2.5.0/rubygems.rb:308:in `activate_bin_path'
/usr/lib/ruby/2.5.0/rubygems.rb:289:in `find_spec_for_exe': can't find gem bundler (>= 0.a) with executable bundler (Gem::GemNotFoundException)
```
    it is because of Gemfile.lock. See https://stackoverflow.com/a/54038218 for
    details. The resolution is by doing the following command:
```
gem install bundler -v 1.16.2
```
1. Pushing the changes on to github will deploy the website.


Use `git tag -a <tag name> -f` to update a tag
