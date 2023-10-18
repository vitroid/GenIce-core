import jinja2 as jj
import toml
import sys

proj = toml.load("pyproject.toml")
t = jj.Template(sys.stdin.read())
print(t.render(**proj))
