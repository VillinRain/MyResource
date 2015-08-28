package parttern.singleton;

public class SingletonDemo {
	/**
	 * 线程安全而且并发度高的单例变种，匿名内部类，匿名了吗？
	 * @author zhaochengyu
	 *http://blog.sina.com.cn/s/blog_75247c770100yxpb.html#bsh-24-289558518
	 */
	private static class  SindletonDemoHandle{
		public final static SingletonDemo instance=new SingletonDemo();
	}
	public static SingletonDemo getInstance(){
		return SindletonDemoHandle.instance;
		
	}
}
