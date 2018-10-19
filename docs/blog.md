---
layout: page
permalink: /blog/
title: "Blog Posts"
---

<ul class="post-list">
  {%- for post in site.posts -%}
  <li>
    {%- assign date_format = site.minima.date_format | default: "%b %-d, %Y" -%}
    <p>
      <span class="post-meta">
        {{ post.date | date: date_format }}
      </span> &mdash;
      <a class="post-link" href="{{ post.url | relative_url }}">
        {{ post.title | escape }}
      </a> 
    </p>
    {%- if site.show_excerpts -%}
      {{ post.excerpt }}
    {%- endif -%}
  </li>
  {%- endfor -%}
</ul>
