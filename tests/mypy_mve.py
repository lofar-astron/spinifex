from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Generic, NamedTuple, Protocol, TypeVar, cast


class FooOptions(NamedTuple):
    foo: str


class BarOptions(NamedTuple):
    bar: str


O = TypeVar("O", FooOptions, BarOptions)  # noqa: E741
O_contra = TypeVar("O_contra", FooOptions, BarOptions, contravariant=True)


class FooBarFunction(Protocol, Generic[O_contra]):
    def __call__(self, options: O_contra | None) -> str: ...


@dataclass
class FooBarFunctions:
    foo: FooBarFunction[FooOptions]
    bar: FooBarFunction[BarOptions]


def foo_func(options: FooOptions | None = None) -> str:
    if options is None:
        return "foo"
    return options.foo


def bar_func(options: BarOptions | None = None) -> str:
    if options is None:
        return "bar"
    return options.bar


foo_bar_functions = FooBarFunctions(
    foo=foo_func,
    bar=bar_func,
)


def parse_options(foo_bar_func: FooBarFunction[O], **kwargs: Any) -> O:
    if foo_bar_func == foo_bar_functions.foo:
        return cast(O, FooOptions(**kwargs))
    if foo_bar_func == foo_bar_functions.bar:
        return cast(O, BarOptions(**kwargs))

    msg = "Unknown function type"
    raise TypeError(msg)


def parse_foo_bar_func(foo_bar_name: str) -> FooBarFunction[O]:
    try:
        return getattr(foo_bar_functions, foo_bar_name)  # type: ignore[no-any-return]
    except AttributeError as e:
        msg = f"Unknown function {foo_bar_name}. Supported functions are {list(foo_bar_functions.__annotations__.keys())}"
        raise TypeError(msg) from e


def baz(
    foo_bar_name: str = "foo",
    **kwargs: Any,
) -> str:
    foo_bar_func = parse_foo_bar_func(foo_bar_name)
    options = parse_options(foo_bar_func, **kwargs)
    return _baz(
        foo_bar_func=foo_bar_func,
        options=options,
    )


def _baz(
    foo_bar_func: FooBarFunction[O],
    options: O | None = None,
) -> str:
    return foo_bar_func(options=options)
